"""Resident Chromnitron model for fast multi-CAP inference (F2 + F3).

Stock load_chromnitron() re-reads the ~340 MB base and re-wraps the LoRA structure once
per CAP. ResidentChromnitron builds base + LoRA ONCE, resident on device, and swaps only
the per-CAP LoRA delta in place. It reuses upstream's helpers (replace_layers_with_lora,
lora_layer_factory, load_state_dict_to_lora) so the per-CAP weight state is byte-identical
to load_chromnitron().

forward_seq() encodes the ESM protein embedding ONCE per CAP and broadcasts it across the
batch (F3). Every norm is GroupNorm/LayerNorm (per-sample; no BatchNorm), so this is
numerically identical to the stock per-batch path. Equivalence: validate_fast_equiv.py.
"""
import os
import torch

from chromnitron_model.chromnitron_models import Chromnitron
from chromnitron_model.load_model import (
    replace_layers_with_lora, lora_layer_factory, load_state_dict_to_lora)


def _load_weights(path, device):
    w = torch.load(path, map_location=device)
    if isinstance(w, dict) and 'model' in w:
        w = w['model']
    return w


class ResidentChromnitron:
    """Build the base+LoRA model once; swap per-CAP LoRA delta in place (F2)."""

    def __init__(self, config, lora_r=4):
        self.config = config
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        mr = config['model_resource']
        self.base_path = f'{mr["root"]}/{mr["base_weights"]}'
        self.lora_dir = f'{mr["root"]}/{mr["per_cap_lora_weights"]}'
        self.finetune_mode = config['inference_config']['inference']['use_finetune']
        assert self.finetune_mode in ['auto', 'enable', 'disable']

        # Cache the base weights in RAM so we never re-read 340 MB from disk per CAP.
        self.base_weights = _load_weights(self.base_path, self.device)

        # Build + wrap ONCE, in the same order load_chromnitron() uses:
        #   fresh model -> load base -> wrap with LoRA -> mark only LoRA trainable.
        import loralib as lora
        self._lora = lora
        model = Chromnitron().to(self.device)
        model.load_state_dict(self.base_weights, strict=False)
        replace_layers_with_lora(model, lora_layer_factory, r=lora_r)
        lora.mark_only_lora_as_trainable(model)
        self.model = model
        self._lora_params = [p for n, p in model.named_parameters() if 'lora_' in n]
        self._cur_cap = None

    def lora_path(self, cap):
        return f'{self.lora_dir}/{cap}.pt'

    def load_lora_state(self, cap, to_host=True):
        """Read a CAP's LoRA .pt (prefetchable to host). Returns None if absent."""
        path = self.lora_path(cap)
        if not os.path.exists(path):
            return None
        return _load_weights(path, 'cpu' if to_host else self.device)

    def set_cap(self, cap, lora_state='__read__'):
        """Make `self.model` equal to stock `load_chromnitron(config, cap)`.

        Reproduces stock's per-CAP result exactly:
          * use_finetune in {auto, enable} and a LoRA file exists -> base + that delta
          * use_finetune == 'enable' but no file                  -> error (as stock)
          * otherwise (no file in 'auto', or 'disable')           -> pure base model

        `lora_state` may be a host-preloaded delta (for prefetch); '__read__' reads now.
        """
        path = self.lora_path(cap)
        weights_exist = os.path.exists(path)
        use_lora = self.finetune_mode in ['auto', 'enable'] and weights_exist
        if self.finetune_mode == 'enable' and not weights_exist:
            raise ValueError(f'Finetuned model for {cap} not found')
        if lora_state == '__read__':
            lora_state = _load_weights(path, self.device) if use_lora else None
        elif not use_lora:
            lora_state = None
        elif lora_state is not None:
            lora_state = {k: v.to(self.device) for k, v in lora_state.items()}

        m = self.model
        # train() lets loralib UNMERGE the previous CAP's delta (W -= s*B@A) using the
        # CURRENT A,B, restoring W to base and clearing the merged flag. The subsequent
        # eval() then re-merges with THIS CAP's delta. Keeping the merge flag in sync is
        # load-bearing for correctness across CAPs.
        m.train()
        with torch.no_grad():
            for p in self._lora_params:
                p.zero_()  # a freshly built model has lora_B=0 -> zero contribution for
                           # any layer this CAP's file does not set; zeroing reproduces that
                           # and prevents a stale delta from leaking between CAPs.
        load_state_dict_to_lora(m, self.base_weights)   # reset conv/norm to base
        if lora_state is not None:
            m.load_state_dict(lora_state, strict=False)  # this CAP's LoRA delta
        m.eval()                                         # merge with this CAP's delta
        self._cur_cap = cap
        return self

    @torch.no_grad()
    def encode_protein(self, esm_embedding):
        """Encode the ESM embedding ONCE per CAP (F3).
        esm_embedding: np.ndarray (prot_len, prot_dim) as stored in <cap>.npz['embedding'].
        Returns (1, hidden, 256), to be broadcast across the batch in forward_seq().
        """
        e = torch.as_tensor(esm_embedding, device=self.device)
        e = e[None].float().transpose(-1, -2)        # (1, prot_dim, prot_len)
        return self.model.prot_encoder(e)            # (1, hidden, 256)

    @torch.no_grad()
    def forward_seq(self, seq, epi, prot_emb):
        """Reproduce Chromnitron.forward for num_targets=1, sample_per_chunk=1, with the
        protein embedding precomputed once and broadcast across the batch.

        seq: (B, 5, 8192) float on device   epi: (B, 1, 8192) float on device
        prot_emb: (1, hidden, 256)          returns preds (B, 8192)
        """
        m = self.model
        seq_embedding = m.encoder((seq, epi))                 # (B, hidden, L)
        B = seq_embedding.size(0)
        H = m.hidden
        L = seq_embedding.size(-1)
        split_emb = torch.zeros(B, H, 1, device=seq_embedding.device, dtype=seq_embedding.dtype)
        prot = prot_emb.expand(B, -1, -1)                     # broadcast single protein emb
        joint = torch.cat([seq_embedding, split_emb, prot], dim=-1)   # (B, H, L+1+256)
        x = joint.transpose(1, 2).contiguous()
        x = m.tf(x)
        x = x.transpose(1, 2).contiguous()
        seq_tf_emb = x[:, :, :L]
        out = m.decoder(seq_tf_emb)
        preds = out[0].view(B, 1, -1)[:, 0, :]                # (B, 8192)
        return preds

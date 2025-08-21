import torch
import torch.nn as nn
import chromnitron_model.chromnitron_blocks as blocks

class Chromnitron(nn.Module):
    def __init__(self, num_genomic_features = 1, hidden = 384, num_attn_blocks = 16, num_of_scale = 4, num_targets = 1, no_confidence_prediction = False, sample_per_chunk = 1, prot_dim = 2560):
        super(Chromnitron, self).__init__()
        print('Initializing Chromnitron')
        encoder_filter_size = 5
        num_blocks = num_of_scale
        self.encoder = blocks.MultiModalEncoder(num_genomic_features, hidden, encoder_filter_size, num_blocks=num_blocks - 1)
        self.tf = blocks.AttentionModule(hidden, num_attn_blocks)
        self.decoder = blocks.Decoder(hidden, hidden, 3, num_blocks=num_blocks, output_channel = num_targets, no_confidence_prediction = no_confidence_prediction)

        self.sample_per_chunk = sample_per_chunk
        self.prot_encoder = blocks.ProteinEncoder(prot_dim, hidden = hidden, filter_size = 5, num_blocks = 3)
        self.hidden = hidden

    def forward(self, seq_feature, prot_feature):
        seq_embedding = self.encoder(seq_feature)
        batch_size = seq_embedding.size(0)
        chunk_size, num_targets, prot_h, prot_len = prot_feature.size()
        assert chunk_size * self.sample_per_chunk == batch_size
        prot_feature = prot_feature.view(chunk_size * num_targets, prot_h, prot_len)
        prot_embedding_chunk = self.prot_encoder(prot_feature)
        prot_embedding_chunk = prot_embedding_chunk.view(chunk_size, num_targets, self.hidden, -1)

        seq_emb_length = seq_embedding.size(-1)

        seq_emb_repeat = seq_embedding.unsqueeze(1).repeat(1, num_targets, 1, 1)
        split_emb = torch.zeros(batch_size, num_targets, self.hidden, 1, device = seq_embedding.device)
        # Repeat the protein embedding within each chunk to save memory
        prot_emb_repeat = prot_embedding_chunk.repeat_interleave(self.sample_per_chunk, 0)
        joint_embedding = torch.cat([seq_emb_repeat, split_emb, prot_emb_repeat], dim = -1)

        joint_batch = joint_embedding.view(batch_size * num_targets, -1, joint_embedding.size(-1))

        x = joint_batch.transpose(1, 2).contiguous()
        x = self.tf(x)
        x = x.transpose(1, 2).contiguous()
        seq_tf_emb = x[:, :, :seq_emb_length]

        out = self.decoder(seq_tf_emb)
        out = [x.view(batch_size, num_targets, -1) for x in out]
        return out
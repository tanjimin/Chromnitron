import torch
import torch.nn as nn
import numpy as np
import copy

class ConvBlock(nn.Module):
    def __init__(self, size, stride = 2, hidden_in = 64, hidden = 64):
        super(ConvBlock, self).__init__()
        pad_len = int(size / 2) 
        self.scale = nn.Sequential(
                        nn.Conv1d(hidden_in, hidden, size, stride, pad_len),
                        nn.GroupNorm(1, hidden),
                        nn.GELU(),
                        )
        self.res = nn.Sequential(
                        nn.Conv1d(hidden, hidden, size, padding = pad_len),
                        nn.GroupNorm(1, hidden),
                        nn.GELU(),
                        nn.Conv1d(hidden, hidden, size, padding = pad_len),
                        nn.GroupNorm(1, hidden),
                        )
        self.relu = nn.ReLU()

    def forward(self, x):
        scaled = self.scale(x)
        identity = scaled
        res_out = self.res(scaled)
        out = self.relu(res_out + identity)
        return out

class ConvTransposeBlock(ConvBlock):
    def __init__(self, size, stride = 2, hidden_in = 64, hidden = 64):
        super().__init__(size, stride, hidden_in, hidden)
        pad_len = size // 2
        self.scale = nn.Sequential(
                        nn.ConvTranspose1d(hidden_in, hidden, size - 1, stride, padding = 0),
                        nn.GroupNorm(1, hidden),
                        nn.GELU(),
                        )
        self.res = nn.Sequential(
                        nn.Conv1d(hidden, hidden, size, padding = pad_len),
                        nn.GroupNorm(1, hidden),
                        nn.GELU(),
                        nn.Conv1d(hidden, hidden, size, padding = pad_len),
                        nn.GroupNorm(1, hidden),
                        )

    def forward(self, x):
        scaled = self.scale(x)
        identity = scaled
        res_out = self.res(scaled)
        out = self.relu(res_out + identity)
        return out

class MultiModalEncoder(nn.Module):
    def __init__(self, num_feat, hidden = 512, filter_size = 9, num_blocks = 2):
        super(MultiModalEncoder, self).__init__()
        hidden_start = hidden // 4
        self.start_seq = ConvBlock(23, 2, hidden_in = 5, hidden = hidden_start)
        self.start_feat = ConvBlock(23, 2, hidden_in = num_feat, hidden = hidden_start)
        hidden_ins = [hidden_start] + [hidden] * (num_blocks - 1)
        hiddens =                     [hidden] * num_blocks
        self.scale_seq = self.get_res_blocks(num_blocks, hidden_ins, hiddens, filter_size)
        self.scale_feat = self.get_res_blocks(num_blocks, hidden_ins, hiddens, filter_size)
        self.conv_end = nn.Conv1d(hidden * 2, hidden, 1)

    def forward(self, x):
        seq, epi = x
        seq = self.scale_seq(self.start_seq(seq))
        epi = self.scale_feat(self.start_feat(epi))

        x = torch.cat([seq, epi], dim = 1)
        out = self.conv_end(x)
        return out

    def get_res_blocks(self, n, his, hs, size):
        blocks = []
        for i, h, hi in zip(range(n), hs, his):
            blocks.append(ConvBlock(size, hidden_in = hi, hidden = h))
        res_blocks = nn.Sequential(*blocks)
        return res_blocks

class ProteinEncoder(nn.Module):
    def __init__(self, num_feat, hidden = 512, filter_size = 9, num_blocks = 4):
        super(ProteinEncoder, self).__init__()
        hidden_start = hidden // 2
        self.start_prot = ConvBlock(23, 2, hidden_in = num_feat, hidden = hidden_start)
        hidden_ins = [hidden_start] + [hidden] * (num_blocks - 1)
        hiddens =                     [hidden] * num_blocks
        self.scale_prot = self.get_res_blocks(num_blocks, hidden_ins, hiddens, filter_size)
        self.conv_end = nn.Conv1d(hidden, hidden, 1)

    def forward(self, x):
        x = self.scale_prot(self.start_prot(x))
        out = self.conv_end(x)
        return out

    def get_res_blocks(self, n, his, hs, size):
        blocks = []
        for i, h, hi in zip(range(n), hs, his):
            blocks.append(ConvBlock(size, hidden_in = hi, hidden = h))
        res_blocks = nn.Sequential(*blocks)
        return res_blocks

class Decoder(nn.Module):
    def __init__(self, in_channel, hidden = 256, filter_size = 3, num_blocks = 2, output_channel = 1, no_confidence_prediction = False):
        super(Decoder, self).__init__()
        hidden_ins = [in_channel] + [hidden] * (num_blocks - 1)
        hiddens =         [hidden] * num_blocks
        self.scale = self.get_res_blocks(num_blocks, hidden_ins, hiddens, filter_size)
        self.conv_end = nn.Conv1d(hidden, output_channel, 1)
        self.confidence_prediction = not no_confidence_prediction
        if self.confidence_prediction:
            self.conv_confidence = nn.Conv1d(hidden, output_channel, 1)

    def forward(self, x):
        x = self.scale(x)
        out = self.conv_end(x)
        if self.confidence_prediction:
            confidence = self.conv_confidence(x)
            return out, confidence
        return out

    def get_res_blocks(self, n, his, hs, size):
        blocks = []
        for i, h, hi in zip(range(n), hs, his):
            blocks.append(ConvTransposeBlock(size, hidden_in = hi, hidden = h))
        res_blocks = nn.Sequential(*blocks)
        return res_blocks

class TransformerLayerPreLN(torch.nn.TransformerEncoderLayer):
    # Pre-LN structure

    def forward(self, src, src_mask = None, src_key_padding_mask = None):
        # MHA section
        src_norm = self.norm1(src)
        src_side, attn_weights = self.self_attn(src_norm, src_norm, src_norm, 
                                    attn_mask=src_mask,
                                    key_padding_mask=src_key_padding_mask)
        src = src + self.dropout1(src_side)

        # MLP section
        src_norm = self.norm2(src)
        src_side = self.linear2(self.dropout(self.activation(self.linear1(src_norm))))
        src = src + self.dropout2(src_side)
        return src

class TransformerEncoder(torch.nn.TransformerEncoder):

    def __init__(self, encoder_layer, num_layers, record_attn = False):
        super(TransformerEncoder, self).__init__(encoder_layer, num_layers)
        self.layers = self._get_clones(encoder_layer, num_layers)
        self.num_layers = num_layers
        #self.norm = norm
        #self.record_attn = record_attn

    def forward(self, src, mask = None, src_key_padding_mask = None):
        output = src
        for mod in self.layers:
            output = mod(output, src_mask=mask, src_key_padding_mask=src_key_padding_mask)
        return output

    def _get_clones(self, module, N):
        return torch.nn.modules.ModuleList([copy.deepcopy(module) for i in range(N)])

class PositionalEncoding(nn.Module):

    def __init__(self, hidden, max_len, dropout = 0.1):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, hidden, 2) * (-np.log(10000.0) / hidden))
        pe = torch.zeros(1, max_len, hidden)
        pe[0, :, 0::2] = torch.sin(position * div_term)
        pe[0, :, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[0, :x.size(1)]
        return self.dropout(x)

class AttentionModule(nn.Module):
    def __init__(self, hidden = 512, layers = 16, max_len = 10000, record_attn = False):
        super(AttentionModule, self).__init__()
        #self.record_attn = record_attn
        self.pos_encoder = PositionalEncoding(hidden, max_len, dropout = 0.1)
        encoder_layers = TransformerLayerPreLN(hidden, 
                                          nhead = 8,
                                          dropout = 0.1,
                                          dim_feedforward = 2048,
                                          batch_first = True)
        self.module = TransformerEncoder(encoder_layers, 
                                         layers, 
                                         record_attn = record_attn)

    def forward(self, x):
        x = self.pos_encoder(x)
        output = self.module(x)
        return output

    def inference(self, x):
        return self.module(x)

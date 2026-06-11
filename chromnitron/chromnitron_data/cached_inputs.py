"""Per-celltype input cache for fast multi-CAP inference (F1).

Stock re-reads and re-transforms the DNA one-hot + ATAC track for every window once per
CAP, but those inputs depend on (window, celltype) only — not the CAP — so a 767-CAP
sweep recomputes them 767x (the dominant CPU cost). CachedCelltypeInputs reads them ONCE
per celltype and reuses them across all CAPs.

Swap-safety: reads go through the same InferenceDataset (Track -> ZarrStorage) as stock,
so a mismatched/swapped DNA or ATAC zarr fails check_data_chr_length() exactly as stock.
The only change is a vectorized one-hot LUT, bit-identical to transforms.to_onehot() (and
additionally accepts upper-case bases — a strict superset).
"""
import numpy as np

from chromnitron_data import transforms

# Vectorized, bit-identical replacement for transforms.to_onehot().
#   a/A -> [1,0,0,0,0]  c/C -> [0,1,0,0,0]  g/G -> [0,0,1,0,0]
#   t/T -> [0,0,0,1,0]  n/N -> [0,0,0,0,1]
_ONEHOT_LUT = np.zeros((128, 5), dtype=np.float32)
_VALID_CP = np.zeros(128, dtype=bool)
for _i, _b in enumerate('acgtn'):
    for _c in (_b, _b.upper()):
        _ONEHOT_LUT[ord(_c), _i] = 1.0
        _VALID_CP[ord(_c)] = True


def to_onehot_fast(seq_u1):
    """seq_u1: numpy '<U1' array (a/c/g/t/n, either case). Returns (L, 5) float32,
    bit-identical to transforms.to_onehot(). Fails loud on any other character so a
    bad base can never become a silent all-zero row."""
    cp = seq_u1.view(np.uint32)
    if cp.size and cp.max() >= 128:
        raise ValueError('unexpected base (code point >= 128) in sequence')
    if cp.size and not _VALID_CP[cp].all():
        bad = sorted({chr(c) for c in cp[~_VALID_CP[cp]].tolist()})
        raise ValueError(f'unexpected base(s) not in ACGTNacgtn: {bad}')
    return _ONEHOT_LUT[cp]


class CachedCelltypeInputs:
    """Read seq one-hot (B,8192,5) + ATAC (B,8192) for every window of `dataset` ONCE.

    `dataset` is an upstream `InferenceDataset` already built for one celltype (its
    `Track`s for seq + ATAC are open and length-checked). Window order matches the
    dataset's region order, i.e. the same order stock produces.
    """

    def __init__(self, dataset, sample_size=8192, n_workers=8):
        region = dataset.region
        n = len(region)
        self.windows = [None] * n
        self.seq = np.zeros((n, sample_size, 5), dtype=np.float32)
        self.atac = np.zeros((n, sample_size), dtype=np.float32)
        seq_track = dataset.data['seq']
        atac_track = dataset.data['input_features']

        def work(i):
            chrom, start_str, end_str, region_id = region[i]
            start, end = int(start_str), int(end_str)
            raw_seq = seq_track.get(chrom, start, end)                 # same read as stock
            self.seq[i] = to_onehot_fast(raw_seq)
            atac = dataset.get_features(atac_track, chrom, start, end)  # same read as stock
            self.atac[i] = transforms.log1p_clip_negative(atac)
            self.windows[i] = (chrom, start, end, region_id)

        # zarr decompression releases the GIL -> threads overlap I/O well.
        if n_workers > 1 and n > 1:
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=n_workers) as ex:
                list(ex.map(work, range(n)))
        else:
            for i in range(n):
                work(i)

    def __len__(self):
        return len(self.windows)

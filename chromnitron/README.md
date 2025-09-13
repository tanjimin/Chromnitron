## Quick Start
* **Google Colab:** Explore Chromnitron's features in an interactive environment [Chromnitron inference colab](https://colab.research.google.com/drive/1NEgKCjJ4ipGMrtq8MuTmWy4LcahGswjv?usp=sharing).

## Local installation guide

1. Clone this package: git clone https://github.com/tanjimin/chromnitron.git
2. Install dependencies:
```bash
conda create -n chromnitron # If you haven't created an environment
conda activate chromnitron
conda install python==3.9
pip install -r requirements_<compute_type>.txt
```

## Download inference data and model weights
1. Download tarballs from dropbox.
```bash
wget -O input_resources.tar \
https://www.dropbox.com/scl/fi/pdtb53vrioxufkyddc5je/input_resources.tar?rlkey=1re0r0k6cisazqyyspcnrhdhg&st=ykb14uen&dl=1

## --- Choose one of the options below --- ##

#This archive includes the base weights and 4 sample LoRA weights, easy for testing <400MB>.
wget -O model_weights_subset.tar \
https://www.dropbox.com/scl/fi/dbi8gc2ej2zi4vktambnn/model_weights_subset.tar?rlkey=gd69we2pq6awjhkfh944iusdz&st=w0i2zgop&dl=1


#This archive includes the base weights and the full set 767 LoRA weights <24GB>.
wget -O model_weights.tar \
https://www.dropbox.com/scl/fi/jlnlydqnp55qzliyj1enb/model_weights.tar?rlkey=j7fypvzyyqhuw36l5gcgmly9s&st=5ixf4p4z&dl=1
```
2. Extract files from archives.
Pick either subset or the full set of model weights: `!tar -xf model_weights_subset.tar` (Subset) or `!tar -xf model_weights.tar` (Fullset)
Sequence, ATAC-seq and CAP embeddings: `!tar -xf input_resources.tar`

1. Move model_weights and input_resources under a common folder following the structure listed below:
```bash
<path-to-chromnitron_resources>
├── input_resources
│   ├── ATAC_seq
│   ├── CAP_embeddings
│   └── DNA_sequence
└── model_weights
    ├── CAPs
    └── chromnitron_base.pt
```

## Edit configuration and prepare for inference
1. Edit inference yaml file configuration at `chromnitron/chromnitron/examples/local_config.yaml`: modify `<path-to-chromnitron_resource>` (Data directory shown above),` <path-to-chromnitron-repository>` (Code repo directory), `<path-to-output-directory>` (Output directory) respectively.

2. Then go to `<path-to-chromnitron-repository>/chromnitron/examples/inputs` and edit the selected genomic regions, cell types (ATAC-seq) and CAPs for inference.

## Run inference:
```bash
python inference.py examples/local_config.yaml
```

# Output file structure
```bash
<path-to-output-directory>
└── HepG2 (Cell type)
    ├── CTCF (CAP1)
    │   ├── output
    │   │   ├── data.npy (Raw inference result)
    │   │   ├── locus.bed (Inference windows [bed])
    │   │   └── locus.csv (Inference windows [csv])
    │   └── processed
    │       ├── data.bigwig (Predicted stored as bigwig)
    │       ├── data.zarr (Prediction stored as zarr)
    │       └── peaks.bed (Peak called on prediction)
    ├── MYC (CAP2)
    │   ├── output
    │   └── processed
    ├── CEBPB (CAP3)
    ...
```

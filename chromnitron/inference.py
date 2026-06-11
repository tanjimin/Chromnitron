# Track inference
import sys
import os
import argparse
import yaml
import numpy as np
import pandas as pd
import torch
from chromnitron_model.load_model import load_chromnitron

def main():
    config_path = sys.argv[1]
    config = load_yaml(config_path)
    loci_info, chrs, celltype_list, cap_list = load_inputs(config) # Load inputs

    # Inference
    if config['inference_config']['inference']['enable'] \
            and config['inference_config']['inference'].get('engine', 'stock') == 'fast':
        # Opt-in throughput engine: removes redundant per-CAP DNA/ATAC/protein loading.
        fast_inference(config, loci_info, chrs, celltype_list, cap_list)
    elif config['inference_config']['inference']['enable']:
        for cap in cap_list:
            print(f'Loading model for {cap}')
            model = load_chromnitron(config, cap)
            device = 'cuda' if torch.cuda.is_available() else 'cpu'
            model.to(device)
            for celltype in celltype_list:
                if verify_prediction_exists(config, celltype, cap): continue
                print(f'Loading data for {celltype}')
                chr_sizes = get_chr_sizes(config, chrs)
                dataloader = load_data(config, celltype, loci_info, cap, chr_sizes)
                print(f'Running inference for {celltype} with {cap}')
                pred_cache, label_df = run_inference(config, model, dataloader, celltype, cap)
                save_prediction(pred_cache, label_df, config, celltype, cap)

    # Post-processing
    if config['inference_config']['post_processing']['enable']:
        import chromnitron_data.postprocessing as postproc
        for cap in cap_list:
            for celltype in celltype_list:
                chr_sizes = get_chr_sizes(config, chrs)
                pred_cache, label_df = load_prediction(config, celltype, cap)
                valid_margin = config['inference_config']['post_processing']['valid_margin']
                data_dict = postproc.pred_to_data_dict(pred_cache, label_df, chr_sizes, valid_margin)
                if config['inference_config']['post_processing']['store_zarr']:
                    postproc.run_store_zarr(config, celltype, cap, data_dict, chr_sizes)
                if config['inference_config']['post_processing']['store_bigwig']:
                    postproc.run_store_bigwig(config, celltype, cap, data_dict, chr_sizes)
                if config['inference_config']['post_processing']['peak_calling']:
                    postproc.run_peak_calling(config, celltype, cap, data_dict)

def load_prediction(config, celltype, cap):
    pred_save_path = f'{config["inference_config"]["output"]["path"]}/{celltype}/{cap}/output/data.npy'
    label_save_path = f'{config["inference_config"]["output"]["path"]}/{celltype}/{cap}/output/locus.csv'
    pred_cache = np.load(pred_save_path)
    label_df = pd.read_csv(label_save_path)
    return pred_cache, label_df

def save_prediction(pred_cache, label_df, config, celltype, cap):
    pred_save_path = f'{config["inference_config"]["output"]["path"]}/{celltype}/{cap}/output/data.npy'
    label_save_path = f'{config["inference_config"]["output"]["path"]}/{celltype}/{cap}/output/locus.csv'
    bed_save_path = f'{config["inference_config"]["output"]["path"]}/{celltype}/{cap}/output/locus.bed'
    os.makedirs(os.path.dirname(pred_save_path), exist_ok=True)
    np.save(pred_save_path, pred_cache)
    label_df.to_csv(label_save_path, index=False)
    label_df[['chr', 'start', 'end', 'region_id']].to_csv(bed_save_path, header=False, index=False, sep='\t')

def load_data(config, celltype, loci_info, cap, chr_sizes):
    input_dict = config['input_resource']
    input_seq_path = os.path.join(input_dict['root'], input_dict['sequence'], f'{config["inference_config"]["input"]["assembly"]}.zarr')
    # Optional atac_root / atac_pattern let the ATAC zarr live on a separate root and/or a
    # nested layout (e.g. '{sample}/s8_normalized/genrich_normalized.zarr'). Both default to
    # the stock behaviour: input_features_path = root/atac/{celltype}.zarr.
    atac_base = input_dict.get('atac_root', input_dict['root'])
    atac_pattern = input_dict.get('atac_pattern')
    atac_rel = atac_pattern.format(sample=celltype) if atac_pattern else f'{celltype}.zarr'
    input_features_path = os.path.join(atac_base, input_dict['atac'], atac_rel)
    assembly = config['inference_config']['input']['assembly']
    esm_feature_path = os.path.join(input_dict['root'], input_dict['cap'], f'{cap}.npz')

    if config['inference_config']['input']['excluded_region_path'] == 'auto':
        excluded_region_path = f"{config['input_resource']['root']}/{config['input_resource']['sequence']}/{assembly}-blacklist.v2.bed"
    else:
        excluded_region_path = config['inference_config']['input']['excluded_region_path']
    if not os.path.exists(excluded_region_path):
        print(f'WARNING: {excluded_region_path} does not exist, using all regions')

    from chromnitron_data.chromnitron_dataset import InferenceDataset
    data = InferenceDataset(loci_info, input_seq_path, input_features_path, esm_feature_path, assembly, chr_sizes, metadata_key = celltype, excluded_region_path = excluded_region_path)

    batch_size = config['inference_config']['inference']['batch_size']
    num_workers = config['inference_config']['inference']['num_workers']
    batch_size = min(batch_size, len(data) // 2 + 1) # Ensure at least 2 batches
    dataloader = torch.utils.data.DataLoader(data, batch_size=batch_size, shuffle=False, num_workers=num_workers)
    return dataloader

def verify_prediction_exists(config, celltype, cap):
    # Define save path
    save_root = config['inference_config']['output']['path']
    save_path = f'{save_root}/{celltype}/{cap}/output'

    # Check if inference data already exists
    anchor_data_path = f'{save_path}/data.npy'
    if os.path.exists(anchor_data_path):
        print(f'Inference data already exists for {save_path}, skipping...')
        return True
    return False

def run_inference(config, model, dataloader, celltype, cap, use_tqdm=True):

    pred_cache = []
    label_cache_dict = {'chr': [],
                        'start': [], 
                        'end': [], 
                        'region_id': [],
                        'celltype' : celltype,
                        'cap' : cap}

    # Use TF32
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True

    # Run inference
    with torch.no_grad():
        if use_tqdm:
            from tqdm import tqdm
            dataloader = tqdm(dataloader)
        for batch in dataloader:
            device = 'cuda' if torch.cuda.is_available() else 'cpu'
            seq, input_features, esm_embeddings, loc_info = batch
            seq = seq.to(device)
            input_features = input_features.to(device)
            esm_embeddings = esm_embeddings.to(device)

            esm_embeddings = esm_embeddings.float().transpose(-1, -2)

            batch_size, mini_bs, seq_len, seq_dim = seq.shape
            seq = seq.view(batch_size * mini_bs, seq_len, seq_dim)
            input_features = input_features.view(batch_size * mini_bs, -1)
            seq = seq.transpose(1, 2).float()
            input_features = input_features.unsqueeze(2).transpose(1, 2).float()

            inputs = (seq, input_features)

            preds, confidence = model(inputs, esm_embeddings)
            preds = preds.detach().cpu().numpy()[:, 0, :]
            pred_cache.append(preds)
            label_cache_dict['start'].extend(loc_info[0].tolist())
            label_cache_dict['end'].extend(loc_info[1].tolist())
            label_cache_dict['chr'].extend(loc_info[2])
            label_cache_dict['region_id'].extend(loc_info[3])

    pred_cache = np.concatenate(pred_cache, axis=0)
    # Exponential transform
    pred_cache = np.exp(pred_cache) - 1
    label_df = pd.DataFrame(label_cache_dict)
    return pred_cache, label_df

def build_fast_label_df(windows, celltype, cap):
    ''' Same columns/order as run_inference's label_df, built from cached windows. '''
    return pd.DataFrame({'chr': [w[0] for w in windows],
                         'start': [w[1] for w in windows],
                         'end': [w[2] for w in windows],
                         'region_id': [w[3] for w in windows],
                         'celltype': celltype,
                         'cap': cap})

def fast_inference(config, loci_info, chrs, celltype_list, cap_list,
                   batch_size=None, cache_workers=8):
    ''' Throughput-optimized inference. Removes three redundancies of the stock loop
    WITHOUT changing the math (see chromnitron_model/resident.py and
    chromnitron_data/cached_inputs.py):
      F1  read + transform seq/ATAC ONCE per celltype, reuse across all CAPs
      F2  build the base + LoRA model ONCE, swap only the per-CAP LoRA delta in place
      F3  encode the protein (ESM) embedding ONCE per CAP, broadcast across the batch
    Reads go through the SAME InferenceDataset/Track/ZarrStorage as stock, so swapping a
    DNA or ATAC zarr is validated identically (chromosome-length check included). '''
    from chromnitron_model.resident import ResidentChromnitron
    from chromnitron_data.cached_inputs import CachedCelltypeInputs

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    on_cuda = device.type == 'cuda'
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True

    bs = batch_size or config['inference_config']['inference']['batch_size']
    input_dict = config['input_resource']
    rm = ResidentChromnitron(config)  # F2: base + LoRA built once, resident on device

    for celltype in celltype_list:
        # Build the dataset ONCE per celltype. This opens AND length-checks the seq/ATAC
        # zarr exactly as stock does; cap_list[0] only supplies an ESM path we never read.
        chr_sizes = get_chr_sizes(config, chrs)
        dataloader = load_data(config, celltype, loci_info, cap_list[0], chr_sizes)
        cache = CachedCelltypeInputs(dataloader.dataset, n_workers=cache_workers)  # F1
        N = len(cache)
        print(f'[fast] {celltype}: cached {N} windows once (reused by all {len(cap_list)} CAPs)')
        if N == 0:
            print(f'[fast] {celltype}: 0 windows after blacklist, skipping')
            continue

        seq_host = torch.from_numpy(cache.seq)
        atac_host = torch.from_numpy(cache.atac)
        if on_cuda:
            seq_host = seq_host.pin_memory()
            atac_host = atac_host.pin_memory()

        for cap in cap_list:
            if verify_prediction_exists(config, celltype, cap):
                continue
            rm.set_cap(cap)  # swap only this CAP's LoRA delta
            esm_path = os.path.join(input_dict['root'], input_dict['cap'], f'{cap}.npz')
            if not os.path.exists(esm_path):
                print(f'[fast] {celltype}/{cap}: missing ESM {esm_path}, skipping')
                continue
            prot_emb = rm.encode_protein(np.load(esm_path)['embedding'])  # F3: once per CAP

            preds_chunks = []
            with torch.no_grad():
                for s0 in range(0, N, bs):
                    seq_b = seq_host[s0:s0 + bs].to(device, non_blocking=on_cuda)
                    atac_b = atac_host[s0:s0 + bs].to(device, non_blocking=on_cuda)
                    seq_in = seq_b.transpose(1, 2).contiguous()   # (b, 5, 8192)
                    epi_in = atac_b.unsqueeze(1)                   # (b, 1, 8192)
                    preds = rm.forward_seq(seq_in, epi_in, prot_emb)
                    preds_chunks.append(preds.float().cpu().numpy())
            pred_cache = np.concatenate(preds_chunks, axis=0)
            pred_cache = np.exp(pred_cache) - 1                    # stock exp transform
            label_df = build_fast_label_df(cache.windows, celltype, cap)
            save_prediction(pred_cache, label_df, config, celltype, cap)
            print(f'[fast] {celltype}/{cap}: done ({N} windows)')

def load_inputs(config):
    loci_info, chrs = read_region_bed(os.path.join(config['inference_config']['input']['root'], config['inference_config']['input']['locus_list_path']))
    celltype_list = read_list(os.path.join(config['inference_config']['input']['root'], config['inference_config']['input']['celltype_list_path']))
    cap_list = read_list(os.path.join(config['inference_config']['input']['root'], config['inference_config']['input']['cap_list_path']))
    return loci_info, chrs, celltype_list, cap_list

def read_region_bed(region_path):
    loci_info = []
    chrs = set()
    with open(region_path, 'r') as file:
        for line in file:
            chrom, start, end = line.strip().split('\t')[:3]
            # Remove any non autosomes
            chr_num = chrom.split('chr')[1]
            if chr_num.isdigit():
                loci_info.append([chrom, int(start), int(end)])
                chrs.add(chrom)
    chrs = list(chrs)
    return loci_info, chrs

def get_chr_sizes(config, chrs):
    chr_sizes = {}
    input_dict = config['input_resource']
    chr_path = os.path.join(input_dict['root'], input_dict['sequence'], f'{config["inference_config"]["input"]["assembly"]}.chrom.sizes')
    with open(chr_path, 'r') as file:
        for line in file:
            chrom, size = line.strip().split('\t')
            chr_sizes[chrom] = int(size)
    # Only includes chromosomes in chrs for easy processing
    chr_sizes = {chr: chr_sizes[chr] for chr in chrs}
    return chr_sizes

def read_list(list_path):
    item_list = []
    with open(list_path, 'r') as file:
        for line in file:
            item_list.append(line.strip())
    return item_list

def load_yaml(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

if __name__ == "__main__":
    main()
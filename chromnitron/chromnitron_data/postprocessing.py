import os
import numpy as np

def pred_to_data_dict(pred_cache, label_df, chr_sizes, valid_margin):
    if valid_margin is None:
        pred_cache = pred_cache[:, valid_margin:-valid_margin]
        label_df['start'] += valid_margin
        label_df['end'] -= valid_margin

    data_dict = {}
    for chr_idx, chr_name in enumerate(chr_sizes):
        chr_track = np.zeros(chr_sizes[chr_name])
        chr_mask = np.zeros(chr_sizes[chr_name])
        for loci_idx, loci_info in label_df.iterrows():
            chr, start, end = loci_info['chr'], loci_info['start'], loci_info['end']
            if chr == chr_name:
                start, end = int(start), int(end)
                # Find out regions where there is overlap
                # Not complete overlap
                if chr_mask[start:end].sum() != 0 and not chr_mask[start:end].sum() == end - start:
                    # Assume overlap is in the front, find the first zero index
                    first_zero_idx = np.where(chr_mask[start:end] == 0)[0][0]
                    # Smoothly merge the data with linear decay
                    old_arr = chr_track[start:start+first_zero_idx]
                    new_arr = pred_cache[loci_idx][:first_zero_idx]
                    # Linear decay
                    old_decayed = old_arr * np.linspace(1, 0, len(old_arr))
                    new_decayed = new_arr * np.linspace(0, 1, len(new_arr))
                    # Merge the data
                    chr_track[start:start+first_zero_idx] = old_decayed + new_decayed
                    # Copy later part of the data
                    chr_track[start+first_zero_idx:end] = pred_cache[loci_idx][first_zero_idx:]
                else:
                    chr_track[start:end] = pred_cache[loci_idx]
                chr_mask[start:end] = 1
        data_dict[chr_name] = chr_track
    return data_dict

def run_peak_calling(config, celltype, cap, data_dict):
    from utils.peak_calling import genome_peak_calling, save_peaks_to_bed
    peak_path = f'{config["inference_config"]["output"]["path"]}/{celltype}/{cap}/processed/peaks.bed'
    if os.path.exists(peak_path):
        print(f'Peak calling for {celltype} with {cap} already exists')
        return
    os.makedirs(os.path.dirname(peak_path), exist_ok=True)
    print(f'Running peak calling for {celltype} with {cap}')
    chrs, starts, ends, peak_values = genome_peak_calling(data_dict)
    save_peaks_to_bed(peak_path, chrs, starts, ends, peak_values)

def run_store_bigwig(config, celltype, cap, data_dict, chr_sizes):
    from utils.io import export_to_bigwig
    bw_path = f'{config["inference_config"]["output"]["path"]}/{celltype}/{cap}/processed/data.bigwig'
    if os.path.exists(bw_path):
        print(f'Bigwig for {celltype} with {cap} already exists')
        return
    os.makedirs(os.path.dirname(bw_path), exist_ok=True)
    print(f'Storing bigwig for {celltype} with {cap}')
    export_to_bigwig(bw_path, chr_sizes, data_dict)

def run_store_zarr(config, celltype, cap, data_dict, chr_sizes):
    from utils.io import export_to_zarr
    zarr_path = f'{config["inference_config"]["output"]["path"]}/{celltype}/{cap}/processed/data.zarr'
    if os.path.exists(zarr_path):
        print(f'Zarr for {celltype} with {cap} already exists')
        return
    os.makedirs(os.path.dirname(zarr_path), exist_ok=True)
    print(f'Storing zarr for {celltype} with {cap}')
    export_to_zarr(zarr_path, chr_sizes, data_dict)

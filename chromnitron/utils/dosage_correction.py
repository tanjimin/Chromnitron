import argparse
import os
import numpy as np
import pandas as pd
import zarr
import numcodecs

def main():
    args = parse_args()
    dosage_correction(args)

def dosage_correction(args):
    data_path = args.data_path
    output_path = args.output_path
    scaler_path = args.scaler_path
    mode = args.mode
    threshold = args.threshold
    # Read scaler_path
    scaler = read_scaler(scaler_path)

    # Convert by scaler
    for cap_idx, row in scaler.iterrows():
        cap_name = row['cap_name']
        scale = float(row['scale'])
        new_dict = {'chrs': {}}
        try:
            data = read_zarr(f'{data_path}/{cap_name}/processed/data.zarr')
            for chr_name in data['chrs'].keys():
                new_dict['chrs'][chr_name] = np.array(data['chrs'][chr_name]) * scale
            os.makedirs(f'{output_path}/{cap_name}/processed/', exist_ok=True)
            save_zarr(f'{output_path}/{cap_name}/processed/data.zarr', new_dict)
        except Exception as e:
            print(f'Error in {cap_name}')
            print(e)

def read_scaler(scaler_path):
    return pd.read_csv(scaler_path, sep=',', header=0)

def read_zarr(zarr_path):
    return zarr.open(zarr_path, mode='r')

def save_zarr(zarr_path, data):
    root = zarr.group(store = zarr_path, overwrite = True) # init zarr group
    zarr_compressor = numcodecs.Blosc(cname = 'zstd', clevel = 3, shuffle = numcodecs.Blosc.SHUFFLE) # Setup compressor
    for chr_name in data['chrs'].keys():
        print('Processing and saving', chr_name)
        chr_arr = data['chrs'][chr_name]
        root.create_dataset(f'chrs/{chr_name}', data = chr_arr, chunks = 1000000, compressor = zarr_compressor)

def parse_args():
    parser = argparse.ArgumentParser()
    # Required arguments:
    # source path, output path, scaler path, mode, and threshold
    parser.add_argument('--data-path', type=str, required=True)
    parser.add_argument('--output-path', type=str, required=True)
    parser.add_argument('--scaler-path', type=str, required=True) # Scaler is in a csv file format of "CAP_name<str>, scale<float>" for each row.
    parser.add_argument('--mode', type=str, required=True, choices=['multiplier', 'filter'])
    parser.add_argument('--threshold', type=float, required=False, default=0.05)
    return parser.parse_args()

if __name__ == '__main__':
    main()
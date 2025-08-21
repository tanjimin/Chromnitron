from multiprocessing import Pool
import sys
import os
import numpy as np
import pandas as pd
import zarr
import numcodecs

def main():
    file_dir = sys.argv[1]
    save_name = sys.argv[2]
    os.makedirs(save_name, exist_ok=True)

    # Read line by line
    lines = []
    current_start = 0
    current_chr = ''
    with open(file_dir, 'r') as f:
        for line in f:
            if not line.startswith('chr'):
                continue
            chr_name, start, end, _ = line.split('\t')
            if chr_name != current_chr: # Switched to new chromosome
                current_chr = chr_name
                current_start = 0
                continue
            if current_start != start: # Fill in gaps
                new_line = f'{chr_name}\t{current_start}\t{start}\t0\n'
                lines.append(new_line)
                lines.append(line)
            else:
                lines.append(line)
            current_start = end

    print('Converting to dataframe...')
    import io
    df = pd.read_csv(io.StringIO('\n'.join(lines)), sep = '\t', header = None)

    root = zarr.group(store = save_name, overwrite = True) # init zarr group
    zarr_compressor = numcodecs.Blosc(cname = 'zstd', clevel = 3, shuffle = numcodecs.Blosc.SHUFFLE) # Setup compressor
    for chr_name in gen_chr_names(df):
        print('Processing and saving', chr_name)
        chr_df = df[df[0] == chr_name]
        if len(chr_df) == 0:
            print(f'No data for {chr_name}')
            continue
        chr_arr = gen_dense_array(chr_df)
        root.create_dataset(f'chrs/{chr_name}', data = chr_arr, chunks = 1000000, compressor = zarr_compressor)

def pool_gen_arr(input_tuple):
    start, end, value = input_tuple
    return np.ones(end - start) * value

def gen_dense_array(chr_df):
    chr_name = chr_df[0].iloc[0]
    array_list = []
    array_list.append(np.zeros(chr_df.iloc[0][1]))
    with Pool(processes=16) as pool:
        new_list = pool.map(pool_gen_arr, chr_df.iloc[:, 1:].values.astype(int))
    array_list.extend(new_list)
    #array_list.append(np.zeros(get_chr_lengths(chr_name) - chr_df.iloc[-1][2]))
    array = np.concatenate(array_list).astype(np.int16)
    return array

def gen_chr_names(df):
    chr_names = df[0].unique()
    # Remove chr with _
    chr_names = [chr_name for chr_name in chr_names if '_' not in chr_name]
    return chr_names

def get_chr_lengths(chr_name):
    raise NotImplementedError('Not implemented')

if __name__ == '__main__':
    main()


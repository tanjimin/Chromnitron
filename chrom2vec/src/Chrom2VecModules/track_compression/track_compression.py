# Compress genomic data with quantization and dynamic binning
import sys
import zarr
import numcodecs
import numpy as np

def main():
    # Load zarr array
    #zarr_path = '/gpfs/data/transcriptronlab/home/jt3545/data/SeqToVec2DB/datalake/zarr_norm/STV88806A90/coverage.zarr'
    #zarr_path = '/gpfs/data/transcriptronlab/home/jt3545/data/SeqToVec2DB/datalake/zarr_norm/STV325D3F16/coverage.zarr'
    zarr_path = sys.argv[1]
    save_zaar_path = sys.argv[2]

    zarr_array = zarr.open(zarr_path, mode='r')

    max_bin_size_power = 4
    random_save_number = np.random.randint(100000)

    chr_dict = {}
    
    for chr_name in zarr_array['chrs'].keys():
        chr_array = zarr_array['chrs'][chr_name]
        background_cutoff = 0.3
        chr_std = qc_std(chr_array, background_cutoff)
        binned_array = dynamic_binning(chr_array, max_bin_size_power)
        chr_dict[chr_name] = [binned_array, {f'qc_std_{background_cutoff}': chr_std}]
    save_zarr(save_zaar_path, chr_dict)

    #save_zarr(f'CTCF_test_{32**max_bin_size_power}_{random_save_number}.zarr', chr_dict)
    #save_bigwig(f'CTCF_test_{32**max_bin_size_power}_{random_save_number}.bw', zarr_array, chr_dict)

def dynamic_binning(signal_array, power_max_bin_size= 10, max_compression_threshold = 0.3, no_compression_threshold = 0.8, edge_buffer = 128, exp_scale = True):
    power_max_bin_size += 1
    array = np.array(signal_array)
    bin_sizes = [32**i for i in range(power_max_bin_size)]
    value_levels = np.linspace(1, 0, len(bin_sizes))
    #log_value_levels = np.log2(value_levels + 1)
    if exp_scale:
        pseudo_num = 2
        value_levels = (np.exp2(value_levels * pseudo_num) - 1) / (2**pseudo_num - 1)
    value_levels = value_levels * no_compression_threshold + (1 - value_levels) * max_compression_threshold
    #print(f'Value levels: {value_levels}')
    bin_index_dict = {}
    for bin_i, (bin_size, exp_value_level) in enumerate(zip(bin_sizes, value_levels)):
        if bin_i == 0:
            bin_index_dict[bin_size] = np.where(array > exp_value_level)[0] 
        elif bin_i == len(bin_sizes) - 1:
            bin_index_dict[bin_size] = np.where(array <= exp_value_level)[0]
        else:
            bin_index_dict[bin_size] = np.where((array > exp_value_level) & (array <= value_levels[bin_i - 1]))[0]
    compressed_array = np.zeros(array.shape)
    for bin_size, bin_index in bin_index_dict.items():
        #print(f'Processing bin size {bin_size}')
        if len(bin_index) == 0:
            print(f'No bin of size {bin_size}')
            continue
        if bin_size == 1:
            compressed_array[bin_index] = array[bin_index]
            continue
        range_starts, range_ends = mask_to_ranges(bin_index)
        bin_starts, bin_ends = split_range_to_bins(range_starts, range_ends, bin_size)
        #print(f'Assigning bin values to {len(bin_starts)} bins')
        #for bin_start, bin_end in zip(tqdm(bin_starts), bin_ends):
        for bin_start, bin_end in zip(bin_starts, bin_ends):
            compressed_array[bin_start:bin_end] = array[bin_start:bin_end].mean()
    return compressed_array

def mask_to_ranges(mask_idx):
    # Turn continuous mask_idx to ranges
    diff = np.diff(mask_idx, prepend=mask_idx[0]-1, append=mask_idx[-1]+1)
    change_points = np.where(diff > 1)[0]
    # insert 0 at the beginning
    change_points = np.insert(change_points, 0, 0)
    # insert last index at the end
    change_points = np.append(change_points, len(mask_idx))
    idx_starts = change_points[:-1]
    idx_ends = change_points[1:] - 1
    range_starts = mask_idx[idx_starts]
    range_ends = mask_idx[idx_ends] + 1
    return range_starts, range_ends

def split_range_to_bins(range_starts, range_ends, bin_size):
    # Split ranges to bins
    bin_starts = []
    bin_ends = []
    #print(f'Processing {len(range_starts)} ranges')
    #for range_start, range_end in zip(tqdm(range_starts), range_ends):
    for range_start, range_end in zip(range_starts, range_ends):
        range_bin_starts = np.arange(range_start, range_end, bin_size).tolist()
        #range_bin_ends = np.append(range_bin_starts[1:], range_end)
        range_bin_ends = range_bin_starts[1:] + [range_end]
        bin_starts.extend(range_bin_starts)
        bin_ends.extend(range_bin_ends)
    bin_starts = np.array(bin_starts)
    bin_ends = np.array(bin_ends)
    return bin_starts, bin_ends

def save_bigwig_chr(bw_name, root, selected_chr_name, signal_array):
    import pyBigWig
    chromsizes = {name: array.shape[0] for name, array in root['chrs'].items()}

    # Define headers
    header = [(k, v) for k, v in chromsizes.items()]

    # Open a new bigWig file for writing
    with pyBigWig.open(bw_name, 'w') as bw:
        # Add the header
        bw.addHeader(header, maxZooms=10)

        # Loop over all chromosomes in the zarr group
        for chr_name, chr_length in header:
            if chr_name != selected_chr_name:
                continue
            print(f'Processing {chr_name}')
            # Load the dense array for this chromosome
            #array = root['chrs'][chr_name][:]
            array = np.array(signal_array)

            diff = np.diff(array, prepend=array[0]-1, append=array[-1]+1)
            change_points = np.where(diff)[0]

            starts = change_points[:-1]
            ends = change_points[1:]
            values = array[starts].astype(float)
            chroms = np.repeat(chr_name, len(starts))

            bw.addEntries(chroms, starts, ends, values=values)


def save_bigwig(bw_name, root, chr_dict):
    import pyBigWig
    chromsizes = {name: array.shape[0] for name, array in root['chrs'].items()}

    # Define headers
    header = [(k, v) for k, v in chromsizes.items()]

    # Open a new bigWig file for writing
    with pyBigWig.open(bw_name, 'w') as bw:
        # Add the header
        bw.addHeader(header, maxZooms=10)

        # Loop over all chromosomes in the zarr group
        for chr_name, chr_length in header:
            print(f'Processing {chr_name}')
            # Load the dense array for this chromosome
            #array = root['chrs'][chr_name][:]
            data, metrics = chr_dict[chr_name]
            array = np.array(data)

            diff = np.diff(array, prepend=array[0]-1, append=array[-1]+1)
            change_points = np.where(diff)[0]

            starts = change_points[:-1]
            ends = change_points[1:]
            values = array[starts].astype(float)
            chroms = np.repeat(chr_name, len(starts))

            bw.addEntries(chroms, starts, ends, values=values)

def save_zarr(zarr_name, chr_dict):
    root = zarr.group(store = zarr_name, overwrite = True) # init zarr group
    zarr_compressor = numcodecs.Blosc(cname = 'zstd', clevel = 3, shuffle = numcodecs.Blosc.SHUFFLE) # Setup compressor
    metrics_list = []
    for chr_name in chr_dict.keys():
        print('Processing and saving', chr_name)
        data, metrics = chr_dict[chr_name]
        chr_arr = np.array(data)
        root.create_dataset(f'chrs/{chr_name}', data = chr_arr, chunks = 1000000, compressor = zarr_compressor)
        # Set attributes for each chromosome
        for key, value in metrics.items():
            root[f'chrs/{chr_name}'].attrs[key] = value
        metrics_list.append(metrics)
    # Set average attributes for the whole group
    for key, value in metrics_list[0].items():
        root.attrs[key] = np.mean([metrics[key] for metrics in metrics_list])

def qc_std(signal_array, cutoff):
    # Calculate variance of each bin
    np.random.seed(0)
    array = np.array(signal_array)
    signal_sample = np.random.choice(array.flatten(), 1000000)
    low_value_index = np.where(signal_sample < cutoff)
    low_value_std = np.std(signal_sample[low_value_index])
    return low_value_std

if __name__ == '__main__':
    main()

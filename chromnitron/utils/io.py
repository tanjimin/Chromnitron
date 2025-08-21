import numpy as np

def export_to_bigwig(bw_name, chr_sizes, data_dict):
    import pyBigWig
    # Define headers
    header = [(k, v) for k, v in chr_sizes.items()]
    # Open a new bigWig file for writing
    with pyBigWig.open(bw_name, 'w') as bw:
        # Add the header
        bw.addHeader(header, maxZooms=10)
        # Loop over all chromosomes in the zarr group
        for chr_name, chr_length in header:
            print(f'Processing {chr_name}')
            # Load the dense array for this chromosome
            array = data_dict[chr_name]
            diff = np.diff(array, prepend=array[0]-1, append=array[-1]+1)
            change_points = np.where(diff)[0]
            starts = change_points[:-1]
            ends = change_points[1:]
            values = array[starts].astype(float)
            chroms = np.repeat(chr_name, len(starts))
            bw.addEntries(chroms, starts, ends, values=values)

def export_to_zarr(zarr_name, chr_sizes, data_dict, chunk_size=1000000):
    import zarr
    import numcodecs
    root = zarr.group(store = zarr_name, overwrite = True) # init zarr group
    zarr_compressor = numcodecs.Blosc(cname = 'zstd', clevel = 3, shuffle = numcodecs.Blosc.SHUFFLE) # Setup compressor
    for chr_name in chr_sizes.keys():
        print('Processing and saving', chr_name)
        chr_arr = data_dict[chr_name]
        root.create_dataset(f'chrs/{chr_name}', data = chr_arr, chunks = chunk_size, compressor = zarr_compressor)

def zarr_to_bigwig(zarr_name, chr_sizes, bigwig_name):
    import zarr
    # Load zarr as data_dict
    data_dict_zarr = zarr.open(zarr_name)
    data_dict = {}
    for chr_name in chr_sizes.keys():
        data_dict[chr_name] = data_dict_zarr[f'chrs/{chr_name}'][:]
    # Export to bigwig
    export_to_bigwig(bigwig_name, chr_sizes, data_dict)

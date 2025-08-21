import sys
import zarr
import pyBigWig
import numpy as np

zarr_name = sys.argv[1]
bw_name = sys.argv[2]

# Open the existing zarr group
root = zarr.open_group(store=zarr_name, mode='r')

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
        array = root['chrs'][chr_name][:]

        diff = np.diff(array, prepend=array[0]-1, append=array[-1]+1)
        change_points = np.where(diff)[0]

        starts = change_points[:-1]
        ends = change_points[1:]
        values = array[starts].astype(float)
        chroms = np.repeat(chr_name, len(starts))

        bw.addEntries(chroms, starts, ends, values=values)

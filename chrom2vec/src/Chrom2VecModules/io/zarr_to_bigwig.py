import argparse
import zarr
import pyBigWig
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Convert a zarr track to bigWig.')
    parser.add_argument('zarr_name')
    parser.add_argument('bw_name')
    parser.add_argument('--clip-negative', action='store_true',
                        help='Clip negative signal values to 0 before writing.')
    parser.add_argument('--log1p', action='store_true',
                        help='Apply log1p to signal values before writing.')
    return parser.parse_args()

def transform_array(array, clip_negative=False, log1p=False):
    # Transformations are opt-in so the default bigWig export remains unchanged.
    if clip_negative:
        array = np.clip(array, 0, np.inf)
    if log1p:
        array = np.log1p(array)
    return array

def main():
    args = parse_args()
    zarr_name = args.zarr_name
    bw_name = args.bw_name

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
            array = transform_array(array, clip_negative=args.clip_negative, log1p=args.log1p)

            diff = np.diff(array, prepend=array[0]-1, append=array[-1]+1)
            change_points = np.where(diff)[0]

            starts = change_points[:-1]
            ends = change_points[1:]
            values = array[starts].astype(float)
            chroms = np.repeat(chr_name, len(starts))

            bw.addEntries(chroms, starts, ends, values=values)

if __name__ == '__main__':
    main()

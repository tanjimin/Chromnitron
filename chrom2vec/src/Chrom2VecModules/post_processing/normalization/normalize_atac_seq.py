import sys
import numpy as np
import pandas as pd
import zarr

def main():
    peak_file = sys.argv[1]
    non_peak_file = sys.argv[2]
    zarr_track_file = sys.argv[3]
    output_file = sys.argv[4]

    # Read peak and non-peak regions (bed file)
    peaks = pd.read_csv(peak_file, sep='\t', header=None)
    non_peaks = pd.read_csv(non_peak_file, sep='\t', header=None)

    # Read zarr track file to numpy array
    zarr_track_dir = zarr.open(zarr_track_file, mode='r')
    zarr_track = {'chrs' : {}}
    for chr_name, track in zarr_track_dir['chrs'].items():
        zarr_track['chrs'][chr_name] = track[:]

    # Calculate the mean signal in non-peak regions
    non_peak_signals = []
    for chr_name, track in zarr_track['chrs'].items():
        non_peaks_chr = non_peaks[non_peaks[0] == chr_name]
        starts = non_peaks_chr.iloc[:, 1]
        ends = non_peaks_chr.iloc[:, 2]
        non_peak_signals_chr = [track[start:end].mean() for start, end in zip(starts, ends)]
        non_peak_signals += non_peak_signals_chr    
    non_peak_mean = np.mean(non_peak_signals)

    # Calculate the 75th percentile signal in peak regions
    peak_signals = []
    for chr_name, track in zarr_track['chrs'].items():
        peaks_chr = peaks[peaks[0] == chr_name]
        starts = peaks_chr.iloc[:, 1]
        ends = peaks_chr.iloc[:, 2]
        peak_signals_chr = [track[start:end].mean() for start, end in zip(starts, ends)]
        peak_signals += peak_signals_chr
    peak_percentile = np.percentile(peak_signals, np.arange(99, 99.9, 0.1)).mean()

    # Initialize the new zarr track similar to the zarr_track_dir
    new_zarr_track_dir = zarr.open_group(output_file, mode='w')
    new_zarr_track = new_zarr_track_dir.create_group('chrs')
    for chr_name, track in zarr_track['chrs'].items():
        new_zarr_track.create_dataset(chr_name, data=track[:], shape=track.shape, dtype='float32')
        
        # Subtract the median signal in non-peak regions
        new_zarr_track[chr_name][:] = new_zarr_track[chr_name][:] - non_peak_mean
        # Divide by the 75th percentile of the signal in peak regions
        new_zarr_track[chr_name][:] = new_zarr_track[chr_name][:] / peak_percentile

if __name__ == '__main__':
    main()

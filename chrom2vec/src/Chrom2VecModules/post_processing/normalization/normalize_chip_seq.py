# The idea is to estimate background noise using non-peak in IP and remove INPUT from IP.

import sys
import numpy as np
import pandas as pd
import zarr

def main():
    peak_file = sys.argv[1]
    non_peak_file = sys.argv[2]
    chip_ip_file = sys.argv[3]
    chip_input_file = sys.argv[4]
    output_file = sys.argv[5]

    # Read peak and non-peak regions (bed file)
    peaks = pd.read_csv(peak_file, sep='\t', header=None)
    non_peaks = pd.read_csv(non_peak_file, sep='\t', header=None)

    # Read zarr track file to numpy array
    chip_ip_track = read_zarr(chip_ip_file)
    chip_input_track = read_zarr(chip_input_file)

    # Calculate the mean signal in non-peak regions
    ip_non_peak_mean = calculate_non_peak_mean(chip_ip_track, non_peaks)
    input_non_peak_mean = calculate_non_peak_mean(chip_input_track, non_peaks)

    # Calculate the correction factor
    correction_factor = ip_non_peak_mean / input_non_peak_mean

    # Apply correction factor to input
    chip_input_track_corrected = {'chrs' : {}}
    for chr_name, track in chip_input_track['chrs'].items():
        chip_input_track_corrected['chrs'][chr_name] = track[:] * correction_factor

    # Correct IP
    chip_ip_track_corrected = {'chrs' : {}}
    for chr_name, track in chip_ip_track['chrs'].items():
        chip_ip_track_corrected['chrs'][chr_name] = track[:] - chip_input_track_corrected['chrs'][chr_name]

    # Adjust IP percentage
    peak_signals = []
    for chr_name, track in chip_ip_track_corrected['chrs'].items():
        peaks_chr = peaks[peaks[0] == chr_name]
        starts = peaks_chr.iloc[:, 1]
        ends = peaks_chr.iloc[:, 2]
        peak_signals_chr = [track[start:end].mean() for start, end in zip(starts, ends)]
        peak_signals += peak_signals_chr

    peak_percentile = mean_value_percentile(np.array(peak_signals), 99, 99.9)
    print(f'Peak percentile: {peak_percentile}')

    for chr_name, track in chip_ip_track_corrected['chrs'].items():
        chip_ip_track_corrected['chrs'][chr_name] = track[:] / peak_percentile

    # Save corrected IP
    save_zarr(chip_ip_track_corrected, output_file)

def mean_value_percentile(peak_signals, low, high):
    low_p, high_p = np.percentile(peak_signals, [low, high])
    p_idx = np.logical_and(peak_signals >= low_p, peak_signals <= high_p)
    if np.sum(p_idx) == 0:
        return np.mean([low_p, high_p])
    return np.mean(peak_signals[p_idx])

def save_zarr(zarr_track, output_file):
    zarr_track_dir = zarr.open_group(output_file, mode='w')
    zarr_track_save = zarr_track_dir.create_group('chrs')
    for chr_name, track in zarr_track['chrs'].items():
        zarr_track_save.create_dataset(chr_name, data=track[:], shape=track.shape, dtype='float32')
        zarr_track_save[chr_name][:] = zarr_track['chrs'][chr_name][:]

def calculate_non_peak_mean(zarr_track, non_peaks):
    # Calculate the mean signal in non-peak regions
    non_peak_signals = []
    for chr_name, track in zarr_track['chrs'].items():
        non_peaks_chr = non_peaks[non_peaks[0] == chr_name]
        starts = non_peaks_chr.iloc[:, 1]
        ends = non_peaks_chr.iloc[:, 2]
        non_peak_signals_chr = [track[start:end].mean() for start, end in zip(starts, ends)]
        non_peak_signals += non_peak_signals_chr    
    non_peak_mean = np.mean(non_peak_signals)
    return non_peak_mean

def read_zarr(zarr_track_file):
    # Read zarr track file to numpy array
    zarr_track_dir = zarr.open(zarr_track_file, mode='r')
    zarr_track = {'chrs' : {}}
    for chr_name, track in zarr_track_dir['chrs'].items():
        zarr_track['chrs'][chr_name] = track[:]
    return zarr_track

if __name__ == '__main__':
    main()


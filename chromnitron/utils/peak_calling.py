import numpy as np
import pandas as pd

def peak_calling(track, threshold):
    '''
    Peak calling function for the inference data. The function will return the peak locations and the peak values.
    '''
    # Find the peaks
    peaks = np.zeros(track.shape)
    peaks[track > threshold] = 1
    diff = np.diff(peaks, prepend=0, append=0)
    change_points = np.where(diff)[0]
    starts = change_points[:-1]
    ends = change_points[1:]
    # Remove peaks that are too small
    peak_values = np.zeros(len(starts))
    for idx, (start, end) in enumerate(zip(starts, ends)):
        peak_values[idx] = np.mean(track[start:end])
    starts = starts[peak_values > threshold]
    ends = ends[peak_values > threshold]
    peak_values = peak_values[peak_values > threshold]
    # Filter out narrow peaks
    narrow_idx = np.where(ends - starts < 100)[0]
    starts = np.delete(starts, narrow_idx)
    ends = np.delete(ends, narrow_idx)
    peak_values = np.delete(peak_values, narrow_idx)
    return starts, ends, peak_values

def genome_peak_calling(data_dict, threshold=0.5):
    '''
    Peak calling function for the genome. The function will return the peak locations and the peak values.
    '''
    # Loop over all chromosomes in the zarr group
    chrs = []
    starts = []
    ends = []
    values = []
    for chr_name in data_dict:
        print(f'Calling peak on {chr_name}')
        # Load the dense array for this chromosome
        track = data_dict[chr_name]
        chr_starts, chr_ends, chr_values = peak_calling(track, threshold)
        chrs.extend([chr_name] * len(chr_starts))
        starts.extend(chr_starts)
        ends.extend(chr_ends)
        values.extend(chr_values)
    return chrs, starts, ends, values

def save_peaks_to_bed(save_path, chrs, starts, ends, values):
    with open(f'{save_path}', 'w') as f:
        for chr, start, end, value in zip(chrs, starts, ends, values):
            f.write(f'{chr}\t{start}\t{end}\t{value}\n')

def load_peaks_from_bed(file_path):
    chrs = []
    starts = []
    ends = []
    values = []
    with open(file_path, 'r') as f:
        for line in f:
            chr, start, end, value = line.strip().split('\t')
            chrs.append(chr)
            starts.append(int(start))
            ends.append(int(end))
            values.append(float(value))
    peak_df = pd.DataFrame({'chr': chrs, 'start': starts, 'end': ends, 'value': values})
    return peak_df

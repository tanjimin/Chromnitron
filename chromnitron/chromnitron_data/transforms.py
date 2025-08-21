import random
import numpy as np

# Onehot
def to_onehot(seq):
    ''' Convert sequence to one-hot encoding '''
    seq_dict = {'a' : [1,0,0,0,0],
                'c' : [0,1,0,0,0],
                'g' : [0,0,1,0,0],
                't' : [0,0,0,1,0],
                'n' : [0,0,0,0,1]}
    onehot_seq = [seq_dict[base] for base in seq]
    return np.array(onehot_seq).astype(np.float32)

def onehot_to_base(onehot_seq):
    ''' Convert one-hot encoding to sequence '''
    seq_dict = {0 : 'a',
                1 : 'c',
                2 : 'g',
                3 : 't',
                4 : 'n'}
    base_seq = [seq_dict[onehot_base] for onehot_base in np.argmax(onehot_seq, axis = 1)]
    return ''.join(base_seq)

# Log(x+1) transform and clip negative values
def log1p_features(input_tracks, output_tracks):
    return log1p_clip_negative(input_tracks), log1p_clip_negative(output_tracks)

def log1p_clip_negative(features):
    ''' Log transform features '''
    if isinstance(features, list):
        return [log1p_clip_negative(f) for f in features]
    else:
        return np.log1p(np.clip(features, 0, np.inf))

# Clip and log transform features
def clip_log_features(input_tracks, output_tracks, clip = 0.001):
    ''' Clip and log transform features '''
    return clipped_log(input_tracks + clip, clip), clipped_log(output_tracks + clip, clip)

def clipped_log(features, minval):
    ''' Clip and log transform features '''
    if isinstance(features, list):
        return [clipped_log(f, minval) for f in features]
    else:
        return np.log(np.clip(features, minval, np.inf))

# Locus subsample and shift
def subsample_locus(window_size, start, end, aug):
    ''' Shift target region '''
    if aug:
        offset = random.choice(range(end - start - window_size))
    else:
        offset = (end - start - window_size) // 2
    return start + offset , start + offset + window_size

# Reverse
def reverse_features(seq, input_tracks, output_tracks, chance = 0.5):
    ''' Reverse augmentation '''
    if random.random() < chance:
        return reverse_complement(seq), reverse_track(input_tracks), reverse_track(output_tracks)
    else:
        return seq, input_tracks, output_tracks

def reverse_track(track):
    ''' Reverse track '''
    if isinstance(track, list):
        return [reverse_track(t) for t in track]
    else:
        return track[::-1].copy()

def reverse_complement(seq):
    ''' Reverse complement onehot vector in form of acgtn '''
    seq = seq[::-1].copy()
    reversed_seq = np.concatenate([seq[:,3:4], seq[:,2:3], seq[:,1:2], seq[:,0:1], seq[:,4:5]], axis = 1)
    return reversed_seq

# Gaussian Noise
def add_gaussian_noise(seq, input_tracks, output_tracks, chance = 0.8):
    ''' Add gaussian noise to input and output tracks '''
    if random.random() < chance:
        return gaussian_noise(seq, std = 0.01), gaussian_noise(input_tracks, std = 0.0001), output_tracks
    else:
        return seq, input_tracks, output_tracks

def gaussian_noise(track, std):
    ''' Add gaussian noise to track'''
    noise = np.random.normal(0, std, track.shape)
    return track + noise

import numpy as np

class Track:
    ''' Genomic track class  '''
    def __init__(self, storage, resolution = 1):
        ''' Initialize track
        storage: storage object
        resolution: storage bin size in base pairs
        '''
        self.storage = storage
        self.resolution = resolution

    def get(self, chrom, start, end):
        ''' Get track data
        chrom: chromosome name
        start: start position
        end: end position
        '''
        storage_start = start // self.resolution
        storage_end = end // self.resolution
        return self.storage.get(chrom, storage_start, storage_end)

    def visualize(self, track_data):
        ''' Visualize track data '''
        raise NotImplementedError

    def save(self, track_data, save_path):
        ''' Save track data '''
        raise NotImplementedError

class AggregatedTrack(Track):
    ''' Aggregated track class '''

    def get(self, chrom, start, end):
        storage_start = start // self.resolution
        storage_end = end // self.resolution
        track_data = [s.get(chrom, storage_start, storage_end) for s in self.storage]
        track_data = np.array(track_data)
        return self.aggregator(track_data)

    def aggregator(self, x):
        raise NotImplementedError

class TrackSum(AggregatedTrack):
    ''' Sum of tracks '''

    def aggregator(self, x):
        return np.sum(x, axis=0)

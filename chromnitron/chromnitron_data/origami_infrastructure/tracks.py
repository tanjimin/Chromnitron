from chromnitron_data.origami_infrastructure.track import Track, TrackSum

class SequenceEmbeddingTrack(Track):

    def visualize(self, track_data):
        track = track_data.mean(axis = 1)

    def save(self, track_data, save_path):
        ''' Save track data '''
        raise NotImplementedError

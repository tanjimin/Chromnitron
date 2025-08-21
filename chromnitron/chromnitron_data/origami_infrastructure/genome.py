
class Genome:
    ''' Abstract geome class for inheritance purpose. '''

    def __init__(self, assembly, chr_sizes):
        self.assembly = assembly
        self.chr_lengths = chr_sizes

    def get_chr_length(self, chr_name):
        ''' Get chromosome length '''
        return self.chr_lengths[chr_name]


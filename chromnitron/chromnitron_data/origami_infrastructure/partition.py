from chromnitron_data.origami_infrastructure.genome import Genome

class Partition(Genome):
    ''' Partition class for a genome. '''

    def __init__(self, loci_info, excluded_loci, assembly, chr_sizes,
                 excluded_chrs=['chrX', 'chrY'],
                 check_length=True, 
                 verbose=False):
        ''' Initialize partition
        loci_info: loci information
        excluded_loci: path to the excluded loci file
        assembly: assembly name
        excluded_chrs: excluded chromosomes
        '''
        super().__init__(assembly, chr_sizes)
        self.loci_info = loci_info
        self.excluded_loci = excluded_loci
        self.verbose = verbose
        self.original_loci = self.load(loci_info) # load partition
        self.loci_on_chrs, self.loci_off_chrs = self.exclude_chrs(self.original_loci) # exclude chromosomes
        self.excluded_loci = self.load_excluded(excluded_loci) # load excluded partition
        self.loci, self.filtered_loci = self.exclude_loci(self.loci_on_chrs, self.excluded_loci) # exclude loci from regions
        if check_length:
            self.check_loci_within_chr()
        else:
            if self.verbose: print('Skip checking loci/chromosome consistency')

    def __len__(self):
        return len(self.loci)

    def __getitem__(self, x):
        return self.loci[x, :]

    def load(self, path):
        ''' Load loci data 
        return: loci, an iterable object
        '''
        raise NotImplementedError

    def load_excluded(self, path):
        ''' Load excluded loci data 
        return: excluded_loci, an iterable object
        '''
        raise NotImplementedError

    def exclude_chrs(self, loci):
        ''' Exclude loci on chromosomes '''
        raise NotImplementedError

    def exclude_loci(self, loci, excluded_loci):
        ''' Exclude loci from the partition '''
        raise NotImplementedError

    def get_loci_chr_location(self, loci_entry):
        ''' Get chromosome name and location from a loci entry '''
        raise NotImplementedError

    def export(self, path):
        ''' Export loci to a bed file '''
        raise NotImplementedError

    def check_loci_within_chr(self):
        ''' Check if the chromosome length in loci is within the chromosome length in the dictionary '''
        for loci_entry in self.loci:
            chr_name, start, end = self.get_loci_chr_location(loci_entry)
            chr_length = self.chr_lengths[chr_name]
            if start < 0 or end > chr_length:
                raise Exception('Loci/chromosome length is not consistent with the chromosome length in the dictionary')
        if self.verbose: print('Loci/chromosome length is consistent.')


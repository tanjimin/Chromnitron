import chromnitron_data.transforms as transforms

from chromnitron_data.origami_infrastructure.tracks import Track
from chromnitron_data.origami_infrastructure.storages import ZarrStorage
from chromnitron_data.origami_infrastructure.partitions import CustomRangeRegion

import numpy as np
from torch.utils.data import Dataset

def get_inference_region(loci_info, assembly, chr_sizes, sample_size, step_size, excluded_region_path, excluded_chrs = ['chrY', 'chrM']):
    region = CustomRangeRegion(sample_size, step_size, loci_info, excluded_region_path, assembly, chr_sizes, excluded_chrs, verbose = False)
    return region

class InferenceDataset(Dataset):
    def __init__(self, loci_info, input_seq_path, input_features_path, esm_feature_path,  # Features
                 assembly, chr_sizes,
                 verbose = False, metadata_key = 'NaN', 
                 sample_size = 8192, step_size = 5120,
                 excluded_region_path = None):
        # Print target features
        if verbose:
            print(f'Loading input seq from {input_seq_path}')
            print(f'Loading input features from {input_features_path}')
        self.metadata_key = metadata_key
        self.input_seq_path = input_seq_path
        self.input_features_path = input_features_path
        self.esm_feature_path = esm_feature_path
        self.verbose = verbose

        self.region = get_inference_region(loci_info, assembly, chr_sizes, sample_size, step_size, excluded_region_path)
        # Initialize data
        self.data = self.load_data(input_seq_path, input_features_path, esm_feature_path, assembly, chr_sizes, verbose)

    def load_data(self, input_seq_path, input_features_path, esm_feature_path, assembly, chr_sizes, verbose):
        ''' Load data from input files '''
        data_dict = {'seq' : None,
                     'input_features' : None,
                     'esm_feature' : None}
        data_dict['seq'] = Track(ZarrStorage(input_seq_path, assembly, chr_sizes))
        data_dict['input_features'] = self.load_storage_with_paths(assembly, 'input_features', input_features_path, chr_sizes)
        data_dict['esm_feature'] = np.load(esm_feature_path)['embedding']
        return data_dict

    def load_storage_with_paths(self, assembly, feature_name, paths, chr_sizes):
        if isinstance(paths, str):
            if paths == '':
                return None
            return Track(ZarrStorage(paths, assembly, chr_sizes))
        elif isinstance(paths, list):
            return [self.load_storage_with_paths(assembly, feature_name, path, chr_sizes) for path in paths]
        else:
            raise ValueError(f'Invalid paths: {paths}')

    def get_features(self, features, chrom, start, end):
        if features is None:
            return np.nan * np.ones((end - start))
        if isinstance(features, Track):
            return features.get(chrom, start, end)
        elif isinstance(features, list):
            return [self.get_features(feature, chrom, start, end) for feature in features]
        else:
            raise ValueError(f'Invalid features: {features}')

    def __len__(self):
        return len(self.region)

    def __getitem__(self, idx):
        # Get sampled region
        chrom, start_str, end_str, region_id = self.region[idx]
        start, end = int(start_str), int(end_str)
        # Get features
        seq = transforms.to_onehot(self.data['seq'].get(chrom, start, end))
        input_features = self.get_features(self.data['input_features'], chrom, start, end)
        esm_feature = self.data['esm_feature']
        # log(1+x) transform features
        input_features = transforms.log1p_clip_negative(input_features)
        # Add zero dimension to esm_feature
        seq = seq[np.newaxis, :, :]
        input_features = input_features[np.newaxis, :]
        esm_feature = esm_feature[np.newaxis, :, :]
        return seq, input_features, esm_feature, (start, end, chrom, region_id, self.metadata_key)

class InferenceSNPDataset(Dataset):

    def __init__(self, dataset, snp_config):
        self.snp_config = snp_config
        self.dataset = dataset
        self.snp_dataset = self.get_snp_dataset(snp_config, dataset)

    def get_snp_dataset(self, snp_config, dataset):
        '''
        loci:
            chr_name: chr5
            location: 1295113
        mutation:
            from: G
            to: A
        '''
        snp_chrom = snp_config['loci']['chr_name']
        snp_location = snp_config['loci']['location']
        wt_base = snp_config['mutation']['from']
        mut_base = snp_config['mutation']['to']

        snp_data_list = []

        for data in dataset:
            seq, input_features, esm_feature, (start, end, chrom, region_id, metadata_key) = data
            # Remove regions that do not contain the SNP
            margin = 2000
            if chrom != snp_chrom or start > snp_location - margin or end < snp_location + margin: continue
            if chrom == snp_chrom and start <= snp_location and end >= snp_location:
                base_seq = transforms.onehot_to_base(seq[0, :, :])
                snp_offset = snp_location - start - 1
                assert base_seq[snp_offset] == wt_base.lower()
                # Covert string to list
                snp_seq_list = list(base_seq)
                snp_seq_list[snp_offset] = mut_base.lower()
                snp_seq = ''.join(snp_seq_list)
                snp_seq = transforms.to_onehot(snp_seq)
                snp_seq = snp_seq[np.newaxis, :, :]
                # Edit region_id
                snp_region_id = f'{region_id}_snp_info_{snp_offset}_{wt_base}_to_{mut_base}'
                snp_data_list.append((snp_seq, input_features, esm_feature, (start, end, chrom, snp_region_id, metadata_key)))
        return snp_data_list

    def __getitem__(self, idx):
        if idx < len(self.dataset):
            return self.dataset[idx]
        else:
            return self.snp_dataset[idx - len(self.dataset)]

    def __len__(self):
        return len(self.dataset) + len(self.snp_dataset)

class InferenceMotifDataset(Dataset):
    def __init__(self, dataset, motif_config = None, mutation_center_list = None):
        self.motif_config = motif_config
        self.mutation_center_list = mutation_center_list
        self.dataset = dataset
        self.motif_dataset = self.get_motif_dataset(motif_config, dataset)

    def __getitem__(self, idx):
        return self.motif_dataset[idx]

    def __len__(self):
        return len(self.motif_dataset)

    def get_motif_dataset(self, motif_config, dataset):
        if self.mutation_center_list is not None:
            mut_center_list = self.mutation_center_list
        else:
            mut_center_list = motif_config['loci']['mutation_center']
        mut_radius = motif_config['loci']['mutation_radius']
        motif_dataset = []
        for data in dataset:
            seq, input_features, esm_feature, (start, end, chrom, region_id, metadata_key) = data
            # Create mutation dataset
            # Check if mutation center is a list
            if not isinstance(mut_center_list, list):
                mut_center = mut_center_list
            else:
                # Select the first mutation center falling into the region
                mut_center = None
                for mut_center_i in mut_center_list:
                    if start <= mut_center_i and end >= mut_center_i:
                        mut_center = mut_center_i
                        break
                if mut_center is None:
                    continue
            mut_start = mut_center - mut_radius
            mut_end = mut_center + mut_radius
            mut_start_offset = mut_start - start
            mut_end_offset = mut_end - start
            seq_list, mut_info_list = self.mutate_seq(seq, mut_start_offset, mut_end_offset)
            for seq, mut_info in zip(seq_list, mut_info_list):
                motif_dataset.append((seq, input_features, esm_feature, (start, end, chrom, region_id + mut_info, metadata_key)))
        return motif_dataset

    def mutate_seq(self, seq, mut_start, mut_end):
        seq_list = []
        seq_list.append(seq)
        mut_info_list = []
        base_seq = transforms.onehot_to_base(seq[0, :, :])
        mut_info_list.append(f'|WT_ref_seq:{base_seq}')
        for i in range(mut_start, mut_end):
            seq_list.extend(self.gen_mutations(seq, i))
            mut_info_list.extend([f'|{i}_mut:A', f'|{i}_mut:C', f'|{i}_mut:G', f'|{i}_mut:T'])
        return seq_list, mut_info_list

    def gen_mutations(self, seq, i):
        seqs = []
        new_seq = seq.copy()
        new_seq[0, i] = [1, 0, 0, 0, 0]
        seqs.append(new_seq)
        new_seq = seq.copy()
        new_seq[0, i] = [0, 1, 0, 0, 0]
        seqs.append(new_seq)
        new_seq = seq.copy()
        new_seq[0, i] = [0, 0, 1, 0, 0]
        seqs.append(new_seq)
        new_seq = seq.copy()
        new_seq[0, i] = [0, 0, 0, 1, 0]
        seqs.append(new_seq)
        return seqs

class InferenceATACPerturbInPlaceDataset(Dataset):

    def __init__(self, dataset, perturb_in_place_config):
        self.perturb_in_place_config = perturb_in_place_config
        self.dataset = dataset
        self.perturb_in_place_dataset = self.get_perturb_in_place_dataset(perturb_in_place_config, dataset)

    def get_perturb_in_place_dataset(self, perturb_in_place_config, dataset):
        '''
        loci:
            chr_name: chr5
            location: 1295113
        mutation:
            from: G
            to: A
        '''
        ptb_chrom = perturb_in_place_config['loci']['chr_name']
        ptb_location = perturb_in_place_config['loci']['location']
        ptb_radius_config = perturb_in_place_config['perturbation_radius']
        ptb_value_range_num = perturb_in_place_config['value_range']
        ptb_value_range = np.arange(ptb_value_range_num[0], ptb_value_range_num[1], ptb_value_range_num[2])

        ptb_data_list = []

        for data in dataset:
            seq, input_features, esm_feature, (start, end, chrom, region_id, metadata_key) = data
            # Remove regions that do not contain the SNP
            margin = 1000
            if chrom != ptb_chrom or start > ptb_location - margin or end < ptb_location + margin: 
                print(f'region {region_id} too far from {ptb_location}')
                continue
            atac_wt = input_features.copy()
            for value in ptb_value_range:
                if isinstance(ptb_radius_config, list):
                    ptb_radius_range = np.arange(ptb_radius_config[0], ptb_radius_config[1], ptb_radius_config[2])
                else:
                    ptb_radius_range = [ptb_radius_config]

                for ptb_radius in ptb_radius_range:
                    ptb_start_relative = ptb_location - start - ptb_radius
                    ptb_end_relative = ptb_location - start + ptb_radius
                    ptb_atac = atac_wt.copy()
                    ptb_atac[:, ptb_start_relative:ptb_end_relative] = atac_wt[:, ptb_start_relative:ptb_end_relative] * value
                    ptb_region_id = f'{region_id}|ptb_atac_info:{value}|peak_center:{ptb_location}|peak_radius:{ptb_radius}'
                    ptb_data_list.append((seq, ptb_atac, esm_feature, (start, end, chrom, ptb_region_id, metadata_key)))
        return ptb_data_list

    def __getitem__(self, idx):
        if idx < len(self.dataset):
            return self.dataset[idx]
        else:
            return self.perturb_in_place_dataset[idx - len(self.dataset)]

    def __len__(self):
        return len(self.dataset) + len(self.perturb_in_place_dataset)
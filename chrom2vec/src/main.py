# A python script that runs the ATAC-seq pipeline on replicates
import subprocess
import os
import time
import yaml
import sys

def main():
    config_path = sys.argv[1]
    config = load_config(config_path)

    # --- Begin execution with time stamp ---
    print(f"Starting pipeline at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")

    fastq_files = config['pipeline_config']['fastq_files']
    save_path = f"{config['pipeline_config']['output_path']}/{config['pipeline_config']['run_name']}"
    module_path = config['pipeline_config']['module_path']
    use_singularity = config['singularity_config']['use_singularity']
    bind_path = config['singularity_config']['bind_path']
    param_tuple = (fastq_files, save_path, module_path, use_singularity, bind_path, config)

    # Softlink, QC and alignment
    softlink_fastq_files(param_tuple)
    run_fastp(param_tuple)
    run_hisat2(param_tuple)

    # Merge and subsample bam files
    merge_bams(param_tuple)
    subsample_bam(param_tuple)

    # Run coverage and post-processing
    run_genrich(param_tuple)
    run_coverage_to_zarr(param_tuple)
    run_normalization(param_tuple) 
    run_zarr_to_bigwig(param_tuple)

# Decorator to calculate runtime for each function and terminal logging
def log_runtime(func):
    def wrapper(*args, **kwargs):
        print(f"\n{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}:\n--- Running {func.__name__} ---")
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        minutes = (end_time - start_time) / 60
        print(f"Function {func.__name__} took {minutes} minutes")
        return result
    return wrapper

def load_config(config_path):
    # Load config file
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

@log_runtime
def softlink_fastq_files(param_tuple):
    fastq_files, save_path, _, _, _, _ = param_tuple
    fastq_path = f'{save_path}/s1_FASTQ'
    # Link fastq files to new directory
    for rep, fastq_reads in fastq_files.items():
        os.makedirs(f"{fastq_path}/{rep}", exist_ok=True)
        for read, fastq_file in fastq_reads.items():
            os.makedirs(f"{fastq_path}/{rep}/{read}", exist_ok=True)
            subprocess.run(["ln", "-s", f"{fastq_file}", f"{fastq_path}/{rep}/{read}/fastq.gz"]) # Soft link to fastq.gz file

@log_runtime
def run_fastp(param_tuple):
    fastq_files, save_path, module_path, use_singularity, bind_path, config = param_tuple
    fastq_path = f'{save_path}/s1_FASTQ'
    fastp_path = f'{save_path}/s2_fastp'

    for rep, fastq_reads in fastq_files.items():
        os.makedirs(f"{fastp_path}/{rep}", exist_ok=True)
        fastp_script_path = f"{module_path}/trimming/fastp_paired_end.sh"
        subprocess.run(["bash", fastp_script_path, use_singularity, bind_path,
                        f"{fastq_path}/{rep}/R1/fastq.gz", f"{fastq_path}/{rep}/R2/fastq.gz",
                        f"{fastp_path}/{rep}/R1.fastq.gz", f"{fastp_path}/{rep}/R2.fastq.gz",
                        f"{fastp_path}/{rep}/fastp.json", f"{fastp_path}/{rep}/fastp.html"])

@log_runtime
def run_hisat2(param_tuple):
    fastq_files, save_path, module_path, use_singularity, bind_path, config = param_tuple
    hisat2_path = f'{save_path}/s3_hisat2'
    fastp_path = f'{save_path}/s2_fastp'

    resources_path = config['pipeline_config']['resources_path']
    hisat2_index = config['resource_config']['hisat2_index']
    assembly_name = config['resource_config']['assembly_name']
    hisat2_index_path = f'{resources_path}/{hisat2_index}/{assembly_name}/genome'
    if not os.path.exists(f'{resources_path}/{hisat2_index}/{assembly_name}'):
        print('Index file does not exist, downloading...')
        download_hisat2_index(param_tuple)
    

    for rep, fastq_reads in fastq_files.items():
        os.makedirs(f"{hisat2_path}/{rep}", exist_ok=True)

        hisat2_script_path = f"{module_path}/alignment/hisat2_paired_end.sh"
        INDEX_PATH = hisat2_index_path
        subprocess.run(["bash", hisat2_script_path, use_singularity, bind_path, INDEX_PATH,
                        f"{fastp_path}/{rep}/R1.fastq.gz", f"{fastp_path}/{rep}/R2.fastq.gz",
                        f"{hisat2_path}/{rep}/hisat2.bam", f"{hisat2_path}/{rep}/summary.txt"])

@log_runtime
def download_hisat2_index(param_tuple):
    _, _, module_path, use_singularity, bind_path, config = param_tuple
    resources_path = config['pipeline_config']['resources_path']
    hisat2_index = config['resource_config']['hisat2_index']
    assembly_name = config['resource_config']['assembly_name']
    hisat2_index_root= f'{resources_path}/{hisat2_index}'
    hisat2_index_download_script_path = f"{module_path}/alignment/hisat2_download_index.sh"
    subprocess.run(["bash", hisat2_index_download_script_path, use_singularity, bind_path, hisat2_index_root, assembly_name])

@log_runtime
def merge_bams(param_tuple):
    fastq_files, save_path, module_path, use_singularity, bind_path, config = param_tuple
    hisat2_path = f'{save_path}/s3_hisat2'
    merged_bam_path = f'{save_path}/s4_merged_bam'
    os.makedirs(f"{merged_bam_path}", exist_ok=True)

    merge_script_path = f"{module_path}/bam/merge_bams.sh"
    inputs = [f"{hisat2_path}/{rep}/hisat2.bam" for rep in fastq_files.keys()]
    subprocess.run(["bash", merge_script_path, use_singularity, bind_path,
                    f"{merged_bam_path}/merged.bam", *inputs])

@log_runtime
def subsample_bam(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    merged_bam_path = f'{save_path}/s4_merged_bam'
    subsampled_bam_path = f'{save_path}/s5_subsampled_bam'
    os.makedirs(f"{subsampled_bam_path}", exist_ok=True)

    subsample_script_path = f"{module_path}/bam/subsample.sh"
    subsample_reads = "40000000"
    subprocess.run(["bash", subsample_script_path, use_singularity, bind_path,
                    "2023", subsample_reads,
                    f"{merged_bam_path}/merged.bam",
                    f"{subsampled_bam_path}/subsampled.bam"])

@log_runtime
def run_genrich(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    subsampled_bam_path = f'{save_path}/s5_subsampled_bam'
    genrich_path = f'{save_path}/s6_genrich'
    os.makedirs(f"{genrich_path}", exist_ok=True)

    genrich_script_path = f"{module_path}/coverage/genrich_coverage.sh"

    resources_path = config['pipeline_config']['resources_path']
    blacklist = config['resource_config']['blacklist']
    blacklist_path = f'{resources_path}/{blacklist}'
    assembly = config['resource_config']['assembly_name']

    os.makedirs(os.path.dirname(blacklist_path), exist_ok=True)
    subprocess.run(["bash", genrich_script_path, use_singularity, bind_path,
                    f"{subsampled_bam_path}/subsampled.bam", blacklist_path, assembly,
                    f"{genrich_path}/genrich.bedgraph"])

@log_runtime
def run_coverage_to_zarr(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    genrich_path = f'{save_path}/s6_genrich'
    zarr_path = f'{save_path}/s7_zarr'
    os.makedirs(f"{zarr_path}", exist_ok=True)

    coverage_script_path = f"{module_path}/io/coverage_to_zarr.sh"
    os.makedirs(f"{zarr_path}", exist_ok=True)
    subprocess.run(["bash", coverage_script_path, use_singularity, bind_path,
                    f"{genrich_path}/genrich.bedgraph",
                    f"{zarr_path}/genrich.zarr"])

@log_runtime
def run_normalization(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    zarr_path = f'{save_path}/s7_zarr'
    normalized_path = f'{save_path}/s8_normalized'
    os.makedirs(f"{normalized_path}", exist_ok=True)


    normalization_script_path = f"{module_path}/post_processing/normalization/normalize_atac_seq.sh"

    resources_path = config['pipeline_config']['resources_path']
    atac_peaks = config['resource_config']['atac_peaks']
    atac_non_peaks = config['resource_config']['atac_non_peaks']
    atac_peaks_path = f'{resources_path}/{atac_peaks}'
    atac_non_peaks_path = f'{resources_path}/{atac_non_peaks}'

    if not os.path.exists(atac_peaks_path):
        print('ATAC peak file does not exist, downloading...')
        peak_download_script_path = f"{module_path}/post_processing/normalization/download_atac_seq_peaks.sh"
        subprocess.run(["bash", peak_download_script_path, use_singularity, bind_path,
                    atac_peaks_path, atac_non_peaks_path])

    subprocess.run(["bash", normalization_script_path, use_singularity, bind_path,
                    atac_peaks_path, atac_non_peaks_path,
                    f"{zarr_path}/genrich.zarr",
                    f"{normalized_path}/genrich_normalized.zarr"])

@log_runtime
def run_zarr_to_bigwig(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    normalized_path = f'{save_path}/s8_normalized'
    bigwig_path = f'{save_path}/s9_bigwig'
    os.makedirs(f"{bigwig_path}", exist_ok=True)

    zarr_to_bigwig_script_path = f"{module_path}/io/zarr_to_bigwig.sh"  

    subprocess.run(["bash", zarr_to_bigwig_script_path, use_singularity, bind_path,
                    f"{normalized_path}/genrich_normalized.zarr",
                    f"{bigwig_path}/genrich_normalized.bw"])

if __name__ == "__main__":
    main()

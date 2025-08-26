from joblib import Parallel, delayed
import os
import pandas as pd
from scipy.stats import gamma
from Bio import SeqIO
import numpy as np
from tqdm import tqdm
import glob
import math
import gzip

def estimate_badread_parameters(read_lengths):
    mean_length = np.mean(read_lengths)
    sd = math.sqrt(np.var(read_lengths))
    return mean_length, sd

def process_single_file(file_path, output_dir):
    result = {}
    if file_path.endswith('.fastq') or file_path.endswith('.fq') or file_path.endswith('.fastq.gz'):
        sample_name = os.path.basename(file_path)
        # Load read lengths from the FASTQ file (handle gzipped files)
        try:
            if file_path.endswith('.gz'):
                with gzip.open(file_path, 'rt') as handle:
                    read_lengths = [len(record.seq) for record in SeqIO.parse(handle, "fastq")]
            else:
                with open(file_path, 'rt') as handle:
                    read_lengths = [len(record.seq) for record in SeqIO.parse(handle, "fastq")]
            # Estimate badread distribution parameters
            mu, sigma = estimate_badread_parameters(read_lengths)
            # Store result
            result = {'Sample': sample_name, 'mean': mu, 'SD': sigma}
        except:
            result = {}
    return result

def process_batch(batch, output_dir):
    batch_results = []
    for file_path in batch:
        result = process_single_file(file_path, output_dir)
        if result:
            batch_results.append(result)
    return batch_results

def process_fastq_files_parallel(sample_list, output_dir, n_jobs=24, batch_size=40):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    # Create batches of files
    batches = [sample_list[i:i + batch_size] for i in range(0, len(sample_list), batch_size)]
    # Process batches in parallel
    results = Parallel(n_jobs=n_jobs)(delayed(process_batch)(batch, output_dir) for batch in tqdm(batches))
    # Flatten list of lists
    results = [result for batch_results in results for result in batch_results if len(result) != 0]
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    return results_df

# Example usage
output_directory = 'simulation_results'  # Replace with the path to your output directory
input_files = glob.glob("Escherichia_coli/nanopore_reads/*")

# Process the FASTQ files in parallel and obtain the results DataFrame
results_df = process_fastq_files_parallel(input_files, output_directory, n_jobs=24, batch_size=100)
results_df.to_csv(os.path.join(output_directory, "parameter_estimates.csv"), index=False)

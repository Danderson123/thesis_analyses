import os
import gzip

def join_fastqs(input_directory, output_file):
    with open(output_file, 'w') as outfile:
        for filename in sorted(os.listdir(input_directory)):
            file_path = os.path.join(input_directory, filename)
            if filename.endswith('.fastq') or filename.endswith('.fq'):
                print(f"Processing {file_path}...")
                with open(file_path, 'r') as infile:
                    for line in infile:
                        outfile.write(line)
            elif filename.endswith('.fastq.gz') or filename.endswith('.fq.gz'):
                print(f"Processing {file_path}...")
                with gzip.open(file_path, 'rt') as infile:
                    for line in infile:
                        outfile.write(line)
    print(f"All FASTQ files have been combined into {output_file}")

# Example usage
if __name__ == "__main__":
    for c in [10, 20, 40, 60, 80, 100, 120, 140]:
        input_dir = f"RW_data/downsampled_reads.{c}x"  # Replace with your input directory path
        output_fastq = f"RW_data/all_samples.{c}x_reads.fastq"  # Replace with your desired output file path
        if not os.path.exists(output_fastq):
            join_fastqs(input_dir, output_fastq)

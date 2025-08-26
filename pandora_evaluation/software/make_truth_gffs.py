import argparse
import glob
import os
import subprocess
parser = argparse.ArgumentParser(description='')
parser.add_argument('--reads-dir', dest='reads_dir', required=True)
parser.add_argument('--reference-fasta', dest='reference_fasta', required=True)
parser.add_argument('--assembly-dir', dest='assembly_dir', required=True)
parser.add_argument('--output-dir', dest='output_dir', required=True)

args = parser.parse_args()

for sample_path in glob.glob(os.path.join(args.reads_dir, "*")):
    sample_name = os.path.basename(sample_path)
    fastq = os.path.join(args.reads_dir, sample_name, sample_name + ".nanopore.fastq.gz")
    assembly = os.path.join(args.assembly_dir, sample_name + ".fa")
    out_gff = os.path.join(args.output_dir, sample_name + ".gff")
    reference_fasta = os.path.join(args.reference_fasta, f"{sample_name}.fasta")
    command = f"python3 software/make_truth_with_minimap.py {assembly} {reference_fasta} {out_gff}"
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    subprocess.run(command, shell=True, check=True)
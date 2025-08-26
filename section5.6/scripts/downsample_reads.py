import os
import subprocess
import glob

all_reads = glob.glob("../amira_paper_data/truth_evaluation/nanopore_reads/SRR32405*")
for c in [10, 20, 40, 60, 80, 100, 120, 140]:
    outdir = f"RW_data/downsampled_reads.{c}x"
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    for reads in all_reads:
        command = f"rasusa reads --seed 2024 --coverage {c} --genome-size 5mb {reads} > {os.path.join(outdir, os.path.basename(reads.replace('.gz', '')))}"
        if not os.path.exists(os.path.join(outdir, os.path.basename(reads.replace('.gz', '')))):
            subprocess.run(command, shell=True, check=True)
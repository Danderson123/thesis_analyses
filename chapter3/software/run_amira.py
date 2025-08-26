import glob
import subprocess
import os
from tqdm import tqdm

fastqs = glob.glob("../pandora_data/*/*.nanopore.fastq.gz")
outdir = "amira_outputs"
if not os.path.exists(outdir):
    os.mkdir(outdir)
for f in tqdm(fastqs):
    sample_output = os.path.join(outdir, os.path.basename(os.path.dirname(f)))
    command = f"singularity run amira.v0.11.0.img amira --no-trim --reads {f} --output {sample_output} --species Escherichia_coli --cores 4 --panRG-path Escherichia.coli.panidx.zip"
    subprocess.run(command, shell=True, check=True)
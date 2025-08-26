import subprocess
import os

if not os.path.exists("metaflye_assemblies"):
    os.mkdir("metaflye_assemblies")
for c in [10, 20, 40, 60, 80, 100, 120, 140]:
    command = f"flye --nano-raw /RW_data/all_samples.{c}x_reads.fastq --meta -t 16 -i 1 -o metaflye_assemblies/RW_data.{c}x_reads"
    if not os.path.exists(f"metaflye_assemblies/RW_data.{c}x_reads"):
        subprocess.run(command, shell=True, check=True)
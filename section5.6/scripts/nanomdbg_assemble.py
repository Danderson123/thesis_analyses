import subprocess
import os

if not os.path.exists("nanomdbg_assemblies"):
    os.mkdir("nanomdbg_assemblies")
for c in [10, 20, 40, 60, 80, 100, 120, 140]:
    command = f"metaMDBG asm --in-ont RW_data/all_samples.{c}x_reads.fastq --threads 16 --out-dir nanomdbg_assemblies/RW_data.{c}x_reads"
    if not os.path.exists(f"nanomdbg_assemblies/RW_data.{c}x_reads"):
        subprocess.run(command, shell=True, check=True)

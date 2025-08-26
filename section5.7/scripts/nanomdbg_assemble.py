import subprocess
import os

if not os.path.exists("ESKAPE_mock_data/nanomdbg_assemblies"):
    os.mkdir("ESKAPE_mock_data/nanomdbg_assemblies")
if not os.path.exists("ESKAPE_mock_data/nanomdbg_amr"):
    os.mkdir("ESKAPE_mock_data/nanomdbg_amr")
for c in ["PF1", "PF2", "PF3", "CC1", "CC2", "CC3", "DM1", "DM2", "DM3", "BS1", "BS2", "BS3"]:
    if not os.path.exists(f"ESKAPE_mock_data/nanomdbg_amr/{c}.tsv"):
        command = f"metaMDBG asm --in-ont ESKAPE_mock_data/ESKAPE_Mock_{c}.fastq.gz --threads 16 --out-dir ESKAPE_mock_data/nanomdbg_assemblies/{c}"
        subprocess.run(command, shell=True, check=True)
        command = f"gunzip ESKAPE_mock_data/nanomdbg_assemblies/{c}/contigs.fasta.gz && amrfinder -d 2024-07-22.1 -n ESKAPE_mock_data/nanomdbg_assemblies/{c}/contigs.fasta -o ESKAPE_mock_data/nanomdbg_amr/{c}.tsv -t 4 --plus"
        subprocess.run(command, shell=True, check=True)

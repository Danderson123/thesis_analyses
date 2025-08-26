import subprocess
import os

if not os.path.exists("ESKAPE_mock_data/metaflye_assemblies"):
    os.mkdir("ESKAPE_mock_data/metaflye_assemblies")
if not os.path.exists("ESKAPE_mock_data/metaflye_amr"):
    os.mkdir("ESKAPE_mock_data/metaflye_amr")
for c in ["PF1", "PF2", "PF3", "CC1", "CC2", "CC3", "DM1", "DM2", "DM3", "BS1", "BS2", "BS3"]:
    command = f"flye --nano-raw ESKAPE_mock_data/ESKAPE_Mock_{c}.fastq.gz --meta -t 16 -i 1 -o ESKAPE_mock_data/metaflye_assemblies/{c}"
    subprocess.run(command, shell=True, check=True)
    command = f"amrfinder -d 2024-07-22.1 -n ESKAPE_mock_data/metaflye_assemblies/{c}/assembly.fasta -o ESKAPE_mock_data/metaflye_amr/{c}.tsv -t 4 --plus"
    subprocess.run(command, shell=True, check=True)

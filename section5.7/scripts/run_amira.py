import os
import subprocess

for sim in ["PF1", "PF2", "PF3", "CC1", "CC2", "CC3", "DM1", "DM2", "DM3", "BS1", "BS2", "BS3"]:
        output = f"ESKAPE_mock_data/amira_output.{sim}"
        threads = 4
        reads = f"ESKAPE_mock_data/ESKAPE_Mock_{sim}.fastq.gz"
        command = f"mkdir -p {output} && singularity run amira.v0.11.0.img amira "
        command += f"--debug --min-relative-depth 0 -n 1 -g 0 --reads {reads} --output {output} --cores {threads} --species metagenome "
        command += f"--amr-fasta AMR_alleles_unified.fa "
        command += "--amr-calls AMR_calls.json "
        command += "--core-genes core_genes.txt "
        command += "--plasmid-genes plasmid_genes.txt "
        command += "--panRG-path ESKAPEES.panidx.zip"
        #if not os.path.exists(os.path.join(output, "amira_results.tsv")):
        subprocess.run(command, shell=True, check=True)
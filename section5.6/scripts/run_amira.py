import os
import subprocess

for sim in ["140x_reads", "120x_reads", "100x_reads", "80x_reads", "60x_reads", "40x_reads", "20x_reads", "10x_reads"]:
        output = f"amira_output/RW_data.{sim}"
        threads = 4
        reads = f"RW_data/all_samples.{sim}.fastq"
        command = f"mkdir -p {output} && singularity run amira.v0.11.0.img amira --debug --meta --reads {reads} --output {output} --species Escherichia_coli --cores {threads} --panRG-path Escherichia.coli.panidx.zip"
        #if not os.path.exists(os.path.join(output, "amira_results.tsv")):
        subprocess.run(command, shell=True, check=True)
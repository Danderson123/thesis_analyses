import os
import subprocess

for sim in ["2", "3", "4", "5", "RW_data"]:
        output = f"/hps/nobackup/iqbal/dander/amira_strain_eval/pandora_output/{sim}"
        pandora_path = "/hps/nobackup/iqbal/dander/Amira_truth_evaluation/software/pandora-linux-precompiled-v0.12.0-alpha.0"
        threads = 32
        panRG = "/hps/nobackup/iqbal/dander/amira_panRG_pipeline/Escherichia_coli_panRG_plasmid_genes/AMR_supplemented_panRG.k15.w5.panidx.zip"
        if sim != "RW_data":
                reads = f"/hps/nobackup/iqbal/dander/amira_strain_eval/strainy_data/ecoli_nano_uni/strain_{sim}/badread_all.fastq.gz"
        else:
                reads = "/hps/nobackup/iqbal/dander/amira_strain_eval/RW_data/all_samples.all_reads.fastq"
        command = f"mkdir -p {output} && {pandora_path} map -t {threads} --debugging-files --min-gene-coverage-proportion 0.5 --max-covg 10000 -o {output} {panRG} {reads}"
        if not os.path.exists(os.path.join(output, "badread_all.fastq.filtered.sam")):
                subprocess.run(command, shell=True, check=True)

"python3 /hps/nobackup/iqbal/dander/Amira_truth_evaluation/scripts/make_truth_with_minimap.py /hps/nobackup/iqbal/dander/amira_strain_eval/metaMDBG_assembly/RW_data.50x_reads/contigs.fasta /hps/nobackup/iqbal/dander/amira_panRG_pipeline/Escherichia_coli_panRG_plasmid_genes/AMR_alleles_unified.fa /hps/nobackup/iqbal/dander/amira_strain_eval/metaMDBG_assembly/RW_data.50x_reads/AMR_genes.json"
import subprocess

command = f"badread simulate --seed 2025 --length 20000,20118 --error_model nanopore2023 --qscore_model nanopore2023 --reference /hps/nobackup/iqbal/dander/amira_metagenome/simulated_metagenome_circular.fna --quantity 1x > /hps/nobackup/iqbal/dander/amira_metagenome/simulated_reads_circular_1.fastq"
subprocess.run(command, shell=True, check=True)
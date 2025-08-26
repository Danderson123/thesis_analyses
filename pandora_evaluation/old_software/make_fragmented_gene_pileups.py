import os
import subprocess

gene = "insB"
outdir = "/hps/nobackup/iqbal/dander/pandora_systematic_evaluation/test"
# remove gaps from the fasta
with open(f"/hps/nobackup/iqbal/dander/pandora_systematic_evaluation/pandora_assessment_with_AMR/0.8_0/panaroo_output/aligned_gene_sequences/{gene}.fasta") as i:
    fasta_content = i.read().split(">")[1:]
cleaned = []
for l in fasta_content:
    header = l.split("\n")[0]
    seq = "".join(l.split("\n")[1:]).replace("-", "")
    cleaned.append(f">{header}\n{seq}")
with open(os.path.join(outdir, f"{gene}.fasta"), "w") as o:
    o.write("\n".join(cleaned))

command = f"minimap2 -x map-ont -t 8 -a --MD --eqx -o {os.path.join(outdir, gene)}_mapped.sam {os.path.join(outdir, f'{gene}.fasta')} {os.path.join(outdir, f'{gene}.fastq.gz')} && samtools sort {os.path.join(outdir, gene)}_mapped.sam > {os.path.join(outdir, gene)}_mapped.bam && samtools index {os.path.join(outdir, gene)}_mapped.bam"
subprocess.run(command, shell=True, check=True)
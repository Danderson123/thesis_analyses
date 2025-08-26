import json
import glob
import os
import random
from tqdm import tqdm

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    comp = ''.join([complement[base] if base in complement else base for base in dna[::-1]])
    return comp

def get_first_allele(fasta_file):
    with open(fasta_file, 'r') as f:
        alignment = f.read().split(">")[1:]
    # Iterate through each line in the FASTA file
    for line in alignment:
        allele_header = line.split("\n")[0]
        sequence = "".join(line.split("\n")[1:]).replace("-", "").upper()
        break
    return os.path.basename(fasta_file).replace(".fasta", "").replace(".aln.fas", ""), allele_header, sequence

def get_alleles(fasta_file):
    with open(fasta_file, 'r') as f:
        alignment = f.read().split(">")[1:]
    # Iterate through each line in the FASTA file
    annotations = {}
    for line in alignment:
        allele_header = line.split("\n")[0]
        sequence = "".join(line.split("\n")[1:]).replace("-", "").upper()
        annotations[allele_header] = sequence
    return os.path.basename(fasta_file).replace(".fasta", "").replace(".aln.fas", ""), annotations

indir = "/hps/nobackup/iqbal/dander/thesis_figures/pandora_systematic_evaluation/pandora_assessment_with_AMR/0.8_0/panaroo_output"
out = "/hps/nobackup/iqbal/dander/thesis_figures/pandora_systematic_evaluation/pandora_assessment_with_AMR/0.8_0/gene_lengths.json"
# list the MSAs
msas = glob.glob(os.path.join(indir, "aligned_gene_sequences", "*.fasta")) + glob.glob(os.path.join(indir, "aligned_gene_sequences", "*.aln.fas"))
random.shuffle(msas)
gene_lengths = {}
for f in tqdm(range(len(msas))):
    # get the header and sequence of the first MSA
    header, annotations = get_alleles(msas[f])
    # iterate through the annotations
    for a in annotations:
        sample = a.split(";")[0].replace("_R_", "")
        # store the gene lengths
        if header not in gene_lengths or len(annotations[a]) < gene_lengths[header]:
            gene_lengths[header] = len(annotations[a])
# write an output file for the gene lengths
with open(out, "w") as o:
    o.write(json.dumps(gene_lengths))

for k in gene_lengths:
    if gene_lengths[k] < 250:
        print(k, gene_lengths[k])
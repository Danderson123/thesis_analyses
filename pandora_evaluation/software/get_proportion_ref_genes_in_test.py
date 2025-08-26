import glob
from tqdm import tqdm
import os

panrg_fastas = glob.glob("pandora_assessment_with_AMR_filtered_with_test/0.8_0/qced_aligned_gene_sequences/*.fasta")
truth_fastas = glob.glob("pandora_assessment_with_AMR_filtered_with_test/0.8_0/panaroo_annotations/GCA*.fasta")

panrg_genes = set()
for f in tqdm(panrg_fastas):
    panrg_genes.add(os.path.basename(f).replace(".fasta", ""))

test_genes = {}
for f in tqdm(truth_fastas):
    with open(f) as i:
        content = i.read().split(">")[1:]
    for allele in content:
        gene = allele.split("\n")[0]
        if gene not in test_genes:
            test_genes[gene] = set()
        test_genes[gene].add(f)

print(len(set(test_genes.keys())), len(panrg_genes))
print(len(set(test_genes.keys()) & panrg_genes))
print(len(set(test_genes.keys()) - panrg_genes))
counts = {}
for g in (set(test_genes.keys()) - panrg_genes):
    if len(test_genes[g]) not in counts:
        counts[len(test_genes[g])] = 0
    counts[len(test_genes[g])] += 1
print(counts)
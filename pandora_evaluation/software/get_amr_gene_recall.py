import glob
import os
from tqdm import tqdm


alns = glob.glob("pandora_assessment_with_AMR/0.95_0.9/panaroo_output/aligned_gene_sequences/*")
amr = set()
for a in tqdm(alns):
    with open(a) as i:
        c = i.read()
    if "AMR_alleles" in c:
        amr.add(os.path.basename(a).replace(".aln.fas", "").replace(".fasta", ""))

read_aln = "pandora_assessment_with_AMR/0.95_0.9/summary_stats.k15.w14_alignments.txt"
with open(read_aln) as i:
    aln_lines = i.read().split(">")[1:]

correct = 0
incorrect = 0
for row in tqdm(aln_lines):
    newline_split = row.split("\n")
    read = newline_split[0]
    truth = newline_split[1].split(",")
    pandora = newline_split[2].split(",")
    for i in range(len(truth)):
        if truth[i][1:] in amr:
            if pandora[i][1:] == truth[i][1:]:
                correct += 1
            else:
                incorrect += 1

print(correct / (correct + incorrect))
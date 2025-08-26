import glob
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import json
from collections import Counter
import numpy as np

# Load read classifications
with open("pandora_assessment_with_AMR_filtered_with_test/0.8_0/reads_longest_contig.txt") as i:
    chromosome = set(i.read().split("\n"))
with open("pandora_assessment_with_AMR_filtered_with_test/0.8_0/reads_other_contigs.txt") as i:
    plasmid = set(i.read().split("\n"))

# Parse annotations
truth_annotations = glob.glob("pandora_assessment_with_AMR_filtered_with_test/0.8_0/processing/*json")

overall_freq = {}

for a in tqdm(truth_annotations):
    with open(a) as i:
        content = json.load(i)
    sample = os.path.basename(a).replace(".json", "")
    for r in content:
        read = r.split("_")[0]
        for g in content[r]:
            gene = g[1:]
            if gene not in overall_freq:
                overall_freq[gene] = set()
            overall_freq[gene].add(sample)

# Count genes per number of unique samples
def frequency_distribution(gene_dict):
    return Counter(len(samples) for samples in gene_dict.values())

overall_counts = frequency_distribution(overall_freq)

# X-axis range
x_vals = sorted(overall_counts)

# Bar heights
overall_heights = [overall_counts.get(x, 0) / len(overall_freq) for x in x_vals]

print(sum(overall_heights))  # Should be ~1

# Bar plot settings
x_indices = np.arange(len(x_vals))

plt.rcParams.update({'font.family': 'sans-serif', 'axes.linewidth': 2})

# Plot only overall
plt.figure(figsize=(10, 8))
plt.bar(x_indices, overall_heights, width=1, color="lightgrey", edgecolor="black")

# Aesthetics
plt.xlabel('Number of samples', fontsize=14)
plt.ylabel('Proportion of genes', fontsize=14)
plt.xticks(x_indices, x_vals, fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0, 0.4])
plt.tight_layout()

# Style
for spine in ['top', 'right']:
    plt.gca().spines[spine].set_visible(False)
plt.tick_params(width=1.2)

# Save
out_path = "/hps/nobackup/iqbal/dander/thesis_figures/intro_pangenome.pdf"
plt.savefig(out_path, dpi=600)
plt.close()

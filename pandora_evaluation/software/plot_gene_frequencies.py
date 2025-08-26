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
c_freq = {}
p_freq = {}

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
            if read in chromosome:
                if gene not in c_freq:
                    c_freq[gene] = set()
                c_freq[gene].add(sample)
            elif read in plasmid:
                if gene not in p_freq:
                    p_freq[gene] = set()
                p_freq[gene].add(sample)
            overall_freq[gene].add(sample)

# Count genes per number of unique samples
def frequency_distribution(gene_dict):
    return Counter(len(samples) for samples in gene_dict.values())

overall_counts = frequency_distribution(overall_freq)
chromosome_counts = frequency_distribution(c_freq)
plasmid_counts = frequency_distribution(p_freq)

# X-axis range
all_x = set(overall_counts) | set(chromosome_counts) | set(plasmid_counts)
x_vals = sorted(all_x)

# Bar heights
overall_heights = [overall_counts.get(x, 0) / len(overall_freq) for x in x_vals]
chromosome_heights = [chromosome_counts.get(x, 0) / len(c_freq)  for x in x_vals]
plasmid_heights = [plasmid_counts.get(x, 0) / len(p_freq) for x in x_vals]

print(sum(overall_heights), sum(chromosome_heights), sum(plasmid_heights))
# Bar plot settings
bar_width = 0.25
x_indices = np.arange(len(x_vals))

# Colorblind-friendly palette
colors = {
    'Overall': '#0072B2',      # Blue
    'Chromosome': '#E69F00',   # Orange
    'Plasmid': '#009E73'       # Green
}

# Plot
plt.figure(figsize=(6, 6))
plt.bar(x_indices - bar_width, overall_heights, width=bar_width, color=colors['Overall'], label='Overall', edgecolor='black')
plt.bar(x_indices, chromosome_heights, width=bar_width, color=colors['Chromosome'], label='Chromosome', edgecolor='black')
plt.bar(x_indices + bar_width, plasmid_heights, width=bar_width, color=colors['Plasmid'], label='Plasmid', edgecolor='black')

# Aesthetics
plt.xlabel('Number of samples', fontsize=14)
plt.ylabel('Proportion of genes', fontsize=14)
plt.xticks(x_indices, x_vals, fontsize=10)
plt.yticks(fontsize=10)
plt.legend(fontsize=12, frameon=False)
plt.tight_layout()

# Style
for spine in ['top', 'right']:
    plt.gca().spines[spine].set_visible(False)
plt.tick_params(width=1.2)

# Save
out_path = "pandora_assessment_with_AMR_filtered_with_test/0.8_0/gene_sample_barplot.png"
plt.savefig(out_path, dpi=600)
plt.close()

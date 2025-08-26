import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import json
from collections import Counter
from matplotlib.patches import Patch
from matplotlib import gridspec
from tqdm import tqdm
import statistics

def extract_metrics(file_path):
    metrics = {
        'recall_overall': None,
        'precision_overall': None,
        'recall_chr': None,
        'precision_chr': None,
        'recall_plasmid': None,
        'precision_plasmid': None,
    }
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("Per gene recall:"):
                metrics['recall_overall'] = float(line.strip().split(":")[1])
            elif line.startswith("Per gene precision:"):
                metrics['precision_overall'] = float(line.strip().split(":")[1])
            elif line.startswith("Chromosome per gene recall:"):
                metrics['recall_chr'] = float(line.strip().split(":")[1])
            elif line.startswith("Chromosome per gene precision:"):
                metrics['precision_chr'] = float(line.strip().split(":")[1])
            elif line.startswith("Plasmid per gene recall:"):
                metrics['recall_plasmid'] = float(line.strip().split(":")[1])
            elif line.startswith("Plasmid per gene precision:"):
                metrics['precision_plasmid'] = float(line.strip().split(":")[1])
    return metrics

def frequency_distribution(gene_dict):
    return Counter(len(samples) for samples in gene_dict.values())

def main(directory):
    with open(os.path.join(directory, "reads_longest_contig.txt")) as f:
        chromosome = set(f.read().splitlines())
    with open(os.path.join(directory, "reads_other_contigs.txt")) as f:
        plasmid = set(f.read().splitlines())

    truth_annotations = glob.glob(os.path.join(directory, "processing", "*json"))
    overall_freq, c_freq, p_freq = {}, {}, {}
    for a in tqdm(truth_annotations):
        with open(a) as f:
            content = json.load(f)
        sample = os.path.basename(a).replace(".json", "")
        for r in content:
            read = r.split("_")[0]
            for g in content[r]:
                gene = g[1:]
                overall_freq.setdefault(gene, set()).add(sample)
                if read in chromosome:
                    c_freq.setdefault(gene, set()).add(sample)
                elif read in plasmid:
                    p_freq.setdefault(gene, set()).add(sample)

    overall_counts = frequency_distribution(overall_freq)
    chromosome_counts = frequency_distribution(c_freq)
    plasmid_counts = frequency_distribution(p_freq)

    x_vals = sorted(set(overall_counts) | set(chromosome_counts) | set(plasmid_counts))
    overall_heights = [overall_counts.get(x, 0) / len(overall_freq) for x in x_vals]
    chromosome_heights = [chromosome_counts.get(x, 0) / len(c_freq) if c_freq else 0 for x in x_vals]
    plasmid_heights = [plasmid_counts.get(x, 0) / len(p_freq) if p_freq else 0 for x in x_vals]
    x_indices = np.arange(len(x_vals))
    bar_width = 0.25

    pattern = os.path.join(directory, "summary_stats.*.k15.w5.txt")
    files = glob.glob(pattern)
    recall_overall, recall_chr, recall_plasmid = [], [], []
    precision_overall, precision_chr, precision_plasmid = [], [], []

    for file in files:
        metrics = extract_metrics(file)
        if all(value is not None for value in metrics.values()):
            recall_overall.append(metrics['recall_overall'])
            recall_chr.append(metrics['recall_chr'])
            recall_plasmid.append(metrics['recall_plasmid'])
            precision_overall.append(metrics['precision_overall'])
            precision_chr.append(metrics['precision_chr'])
            precision_plasmid.append(metrics['precision_plasmid'])

    data = [
        recall_overall, recall_chr, recall_plasmid,
        precision_overall, precision_chr, precision_plasmid
    ]
    print(statistics.mean(recall_overall), statistics.mean(recall_chr), statistics.mean(recall_plasmid))
    print(statistics.mean(precision_overall), statistics.mean(precision_chr), statistics.mean(precision_plasmid))

    positions = [0.8, 1.0, 1.2, 1.8, 2.0, 2.2]

    colors = {
        'Overall': '#0072B2',
        'Chromosome': '#E69F00',
        'Plasmid': '#009E73'
    }
    group_colors = [
        colors['Overall'], colors['Chromosome'], colors['Plasmid'],
        colors['Overall'], colors['Chromosome'], colors['Plasmid']
    ]
    legend_patches = [
        Patch(facecolor=colors['Overall'], edgecolor='black', label='All genes'),
        Patch(facecolor=colors['Chromosome'], edgecolor='black', label='Chromosome genes'),
        Patch(facecolor=colors['Plasmid'], edgecolor='black', label='Plasmid genes')
    ]

    # Updated layout: 2 rows, 2 columns (second column for legend)
    fig = plt.figure(figsize=(10, 12))
    spec = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[4, 1], height_ratios=[1, 1])

    # Barplot (top left)
    ax1 = fig.add_subplot(spec[0, 0])
    ax1.bar(x_indices - bar_width, overall_heights, width=bar_width, color=colors['Overall'], edgecolor='black')
    ax1.bar(x_indices, chromosome_heights, width=bar_width, color=colors['Chromosome'], edgecolor='black')
    ax1.bar(x_indices + bar_width, plasmid_heights, width=bar_width, color=colors['Plasmid'], edgecolor='black')
    ax1.set_xlabel('Number of samples', fontsize=12)
    ax1.set_ylabel('Proportion of genes', fontsize=12)
    ax1.set_xticks(x_indices)
    ax1.set_xticklabels(x_vals, fontsize=10)
    ax1.tick_params(width=1.5)
    for spine in ['left', 'bottom']:
        ax1.spines[spine].set_linewidth(1.2)
    for spine in ['top', 'right']:
        ax1.spines[spine].set_visible(False)

    # Boxplot (bottom left)
    ax2 = fig.add_subplot(spec[1, 0])
    boxplots = ax2.boxplot(data,
                           positions=positions,
                           widths=0.18,
                           patch_artist=True,
                           showfliers=False,
                           boxprops=dict(linewidth=1.2),
                           medianprops=dict(color='black', linewidth=1.5),
                           whiskerprops=dict(color='black', linewidth=1.2),
                           capprops=dict(color='black', linewidth=1.2))
    for patch, color in zip(boxplots['boxes'], group_colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
    for i, values in enumerate(data):
        x = np.random.normal(loc=positions[i], scale=0.015, size=len(values))
        ax2.scatter(x, values, color='black', s=15, alpha=0.7, zorder=3)
    ax2.set_xticks([1.0, 2.0])
    ax2.set_xticklabels(['Recall', 'Precision'], fontsize=10)
    ax2.set_ylabel('Value', fontsize=12)
    ax2.set_xlim(0.5, 2.5)
    ax2.set_ylim(0.5, 1.0)
    ax2.tick_params(width=1.5)
    for spine in ['left', 'bottom']:
        ax2.spines[spine].set_linewidth(1.2)
    for spine in ['top', 'right']:
        ax2.spines[spine].set_visible(False)

    # Legend (right middle)
    ax_legend = fig.add_subplot(spec[:, 1])
    ax_legend.axis('off')
    ax_legend.legend(
        handles=legend_patches,
        loc='center',
        fontsize=12,
        frameon=False
    )

    plt.tight_layout()
    out_path = os.path.join(directory, "test_results.png")
    plt.savefig(out_path, dpi=600)
    plt.savefig(out_path.replace(".png", ".pdf"))
    plt.close()

if __name__ == "__main__":
    main("pandora_assessment_with_AMR_filtered_with_test/0.8_0")

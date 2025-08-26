import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

def extract_metrics(file_path):
    """Extract recall and precision metrics from a stats file."""
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

def main(directory):
    pattern = os.path.join(directory, "summary_stats.*.k15.w5.txt")
    files = glob.glob(pattern)

    recall_overall = []
    recall_chr = []
    recall_plasmid = []
    precision_overall = []
    precision_chr = []
    precision_plasmid = []

    for file in files:
        metrics = extract_metrics(file)
        if all(value is not None for value in metrics.values()):
            recall_overall.append(metrics['recall_overall'])
            recall_chr.append(metrics['recall_chr'])
            recall_plasmid.append(metrics['recall_plasmid'])
            precision_overall.append(metrics['precision_overall'])
            precision_chr.append(metrics['precision_chr'])
            precision_plasmid.append(metrics['precision_plasmid'])

    # Data and positions
    data = [
        recall_overall, recall_chr, recall_plasmid,
        precision_overall, precision_chr, precision_plasmid
    ]
    positions = [0.8, 1.0, 1.2, 1.8, 2.0, 2.2]

    # Colorblind-friendly palette
    colors = {
        'Overall': '#0072B2',     # Blue
        'Chromosome': '#E69F00', # Orange
        'Plasmid': '#009E73'     # Green
    }

    group_colors = [
        colors['Overall'], colors['Chromosome'], colors['Plasmid'],
        colors['Overall'], colors['Chromosome'], colors['Plasmid']
    ]

    # Plotting
    plt.figure(figsize=(6, 6))
    boxplots = plt.boxplot(data,
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

    # Jittered points
    for i, values in enumerate(data):
        x = np.random.normal(loc=positions[i], scale=0.015, size=len(values))
        plt.scatter(x, values, color='black', s=15, alpha=0.7, zorder=3)

    # X-axis: single labels
    plt.xticks([1.0, 2.0], ['Recall', 'Precision'], fontsize=10)

    # Legend
    legend_patches = [
        Patch(facecolor=colors['Overall'], edgecolor='black', label='Overall'),
        Patch(facecolor=colors['Chromosome'], edgecolor='black', label='Chromosome'),
        Patch(facecolor=colors['Plasmid'], edgecolor='black', label='Plasmid')
    ]
    plt.legend(handles=legend_patches, fontsize=11, frameon=False, loc='lower left')

    # Style
    plt.yticks(fontsize=10)
    plt.ylabel('Score', fontsize=12)
    plt.xlim(0.5, 2.5)
    plt.ylim(0.5, 1.0)
    plt.tick_params(width=1.2)
    for spine in ['top', 'right']:
        plt.gca().spines[spine].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(directory, "grouped_boxplot.png"), dpi=600)
    plt.close()

if __name__ == "__main__":
    main("/hps/nobackup/iqbal/dander/thesis_figures/pandora_systematic_evaluation/pandora_assessment_with_AMR_filtered_with_test/0.8_0")

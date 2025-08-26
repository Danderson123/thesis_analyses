import matplotlib.pyplot as plt
import os
import glob
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.cm as cm
import statistics

# Set Nature-style but sans-serif aesthetics
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "black",
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "legend.frameon": False
})

# Helper function to gather file paths
def get_files(path):
    return [
        f for f in glob.glob(path)
        if not ("_alignments.txt" in f or ".fp_genes.txt" in f or "H131800734" in f)
    ]

# Collect files for each stage
no_correction = get_files("assessment_results/summary_stats_no_correction/*.txt")
pre_correction = get_files("assessment_results/summary_stats_pre_correction/*.txt")
mid_correction = get_files("assessment_results/summary_stats_mid_correction/*.txt")
post_correction = get_files("assessment_results/summary_stats_post_correction/*.txt")

# Function to extract precision and recall
def extract_stats(file_list, stage):
    records = []
    recalls = []
    precisions = []
    for f in file_list:
        with open(f) as i:
            lines = i.read().splitlines()
        sample = os.path.basename(f).replace(".txt", "")
        recall = next((float(line.split(": ")[1]) for line in lines if "Per gene recall" in line), None)
        precision = next((float(line.split(": ")[1]) for line in lines if "Per gene precision" in line), None)
        records.append({"sample": sample, "recall": recall, "precision": precision, "stage": stage})
        recalls.append(recall)
        precisions.append(precision)
    #print(f"{stage} {statistics.mean(recalls)} {statistics.mean(precisions)}")
    return records

# Combine data
stage_order = ["no correction", "partial gene filtering", "node filtering", "iterative correction"]
all_data = (
    extract_stats(no_correction, stage_order[0]) +
    extract_stats(pre_correction, stage_order[1]) +
    extract_stats(mid_correction, stage_order[2]) +
    extract_stats(post_correction, stage_order[3])
)

df = pd.DataFrame(all_data)

# Ensure correct data types
df["recall"] = pd.to_numeric(df["recall"], errors='coerce')
df["precision"] = pd.to_numeric(df["precision"], errors='coerce')
df["stage"] = pd.Categorical(df["stage"], categories=stage_order, ordered=True)

# Assign a unique color to each sample using a colormap
samples = sorted(df["sample"].unique())
color_map = cm.get_cmap("tab20", len(samples))
sample_colors = {sample: color_map(i) for i, sample in enumerate(samples)}

# Prepare figure
fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

def stylize_boxplot(ax):
    # Force black edges and white fill on boxes
    for patch in ax.patches:
        patch.set_edgecolor('black')
        patch.set_facecolor('white')
        patch.set_linewidth(1)

    # Set black for all lines: medians, whiskers, caps
    for line in ax.lines:
        line.set_color('black')
        line.set_linewidth(1)

# For legend handles
legend_handles = []

# --- Recall plot ---
sns.boxplot(data=df, x="stage", y="recall", ax=axes[0], showfliers=False, width=0.6)
stylize_boxplot(axes[0])

for i, (sample, subdf) in enumerate(df.groupby("sample")):
    line, = axes[0].plot(subdf["stage"].astype(str).to_numpy(), subdf["recall"].to_numpy(),
                         color=sample_colors[sample], marker="o", linewidth=0.8, markersize=3)
    # Only add to legend once
    legend_handles.append((line, sample))

axes[0].set_ylabel("Recall", fontsize=10)
axes[0].set_ylim(0.97, 1.0)
axes[0].yaxis.set_major_locator(MultipleLocator(0.01))
axes[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# --- Precision plot ---
sns.boxplot(data=df, x="stage", y="precision", ax=axes[1], showfliers=False, width=0.6)
stylize_boxplot(axes[1])

for sample, subdf in df.groupby("sample"):
    axes[1].plot(subdf["stage"].astype(str).to_numpy(), subdf["precision"].to_numpy(),
                 color=sample_colors[sample], marker="o", linewidth=0.8, markersize=3)

axes[1].set_ylabel("Precision", fontsize=10)
axes[1].set_xlabel("Correction stage", fontsize=10)
axes[1].set_ylim(0.97, 1.0)
axes[1].yaxis.set_major_locator(MultipleLocator(0.01))
axes[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Add legend to center right
handles, labels = zip(*legend_handles)
fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.4, 0.5))

for ax in axes:
    ax.tick_params(axis='x', labelrotation=25)
    ax.tick_params(axis='both', length=3, width=0.8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

# Optional: Check drops
recall_pairs = {}
precision_pairs = {}
for index, row in df.iterrows():
    if row["stage"] not in ["node filtering", "iterative correction"]:
        continue
    if row["sample"] not in recall_pairs:
        recall_pairs[row["sample"]] = {}
    recall_pairs[row["sample"]][row["stage"]] = row["recall"]
    if row["sample"] not in precision_pairs:
        precision_pairs[row["sample"]] = {}
    precision_pairs[row["sample"]][row["stage"]] = row["precision"]

for s in precision_pairs:
    if precision_pairs[s]["iterative correction"] < precision_pairs[s]["node filtering"]:
        print(s)

plt.tight_layout()
output_path = "assessment_results/correction_trajectories.png"
plt.savefig(output_path, dpi=600, bbox_inches='tight')
plt.savefig(output_path.replace(".png", ".pdf"), dpi=600, bbox_inches='tight')
plt.close()
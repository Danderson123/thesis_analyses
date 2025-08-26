import glob
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import gzip
import json
import seaborn as sns
from matplotlib.lines import Line2D

def parse_fastq_lines(fh):
    # Initialize a counter to keep track of the current line number
    line_number = 0
    # Iterate over the lines in the file
    for line in fh:
        # Increment the line number
        line_number += 1
        # If the line number is divisible by 4, it's a sequence identifier line
        if line_number % 4 == 1:
            # Extract the identifier from the line
            identifier = line.split(" ")[0][1:]
        # If the line number is divisible by 4, it's a sequence line
        elif line_number % 4 == 2:
            sequence = line.strip()
        elif line_number % 4 == 0:
            # Yield the identifier, sequence and quality
            yield identifier, sequence, line.strip()


def parse_fastq(fastq_file):
    # Initialize an empty dictionary to store the results
    results = {}
    # Open the fastq file
    if ".gz" in fastq_file:
        try:
            with gzip.open(fastq_file, "rt") as fh:
                # Iterate over the lines in the file
                for identifier, sequence, quality in parse_fastq_lines(fh):
                    # Add the identifier and sequence to the results dictionary
                    results[identifier.replace("\n", "")] = {
                        "sequence": sequence,
                        "quality": quality,
                    }
            return results
        except OSError:
            pass
    with open(fastq_file, "r") as fh:
        # Iterate over the lines in the file
        for identifier, sequence, quality in parse_fastq_lines(fh):
            # Add the identifier and sequence to the results dictionary
            results[identifier.replace("\n", "")] = {"sequence": sequence, "quality": quality}
    # Return the dictionary of results
    return results

mean_read_lengths = ["5", "10", "20", "40"]
mean_read_depths = ["5", "10", "20", "40", "80"]

if not os.path.exists("simulation_results/read_lengths_aggregated.json"):
    all_read_lengths = {}
    for m in tqdm(mean_read_lengths):
        all_read_lengths[m] = {}
        for d in mean_read_depths:
            sim_directory = f"simulations_{m}kb_r10"
            fastqs = glob.glob(os.path.join(sim_directory, f"reads_{d}", "*.fastq"))
            per_read_lengths = []
            for f in tqdm(fastqs):
                reads = parse_fastq(f)
                for r in reads:
                    per_read_lengths.append(len(reads[r]["sequence"]))
            all_read_lengths[m][d] = per_read_lengths
    with open("simulation_results/read_lengths_aggregated.json", "w") as o:
        o.write(json.dumps(all_read_lengths))
with open("simulation_results/read_lengths_aggregated.json") as i:
    all_read_lengths = json.load(i)

# Font and layout settings
plt.rcParams.update({
    'font.family': 'sans-serif',  # Nature uses sans-serif (e.g., Arial/Helvetica)
    'font.size': 12,
    'axes.linewidth': 0.5,
    'axes.color': "black",
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'xtick.direction': 'out',
    'ytick.direction': 'out'
})

# Create subplots
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(7, 6.5), constrained_layout=True, sharex=True, sharey=True)
axes = axes.flatten()

# Colorblind-friendly palette
palette = sns.color_palette("colorblind", len(mean_read_depths))[::-1]

# Plot histograms
for idx, m in enumerate(mean_read_lengths):
    ax = axes[idx]
    for depth_idx, d in enumerate(reversed(mean_read_depths)):
        ax.hist(
            all_read_lengths[m][d],
            bins=50,
            label=f"{d}×",
            histtype='step',  # Outline only for clarity
            linewidth=2,
            log=True,
            color=palette[depth_idx],
            alpha=1
        )
        print(m, d, max(all_read_lengths[m][d]))
    # Subplot title
    ax.set_title(f"$\mu$ = {m} kb", fontsize=12)

    # Axis labels
    if idx in [2, 3]:
        ax.set_xlabel("Read length (bp)", fontsize=12)
    if idx in [0, 2]:
        ax.set_ylabel("Count (log scale)", fontsize=12)

    # Axis style
    #ax.grid(True, linestyle="-", alpha=0.3, linewidth=0.8)
    ax.set_xlim([0, 400000])
    ax.set_ylim([1, 1e5])
    ax.set_xticks([0, 200000, 400000])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color("black")
    ax.spines['left'].set_color("black")
    #ax.spines['left'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)

# Create custom line handles for the legend
# Reverse both colors and depths for legend consistency
reversed_depths = list(reversed(mean_read_depths))
legend_lines = [
    Line2D([0], [0], color=palette[i], linewidth=2, linestyle='-') 
    for i in range(len(reversed_depths))
]
legend_labels = [f"{d}×" for d in reversed_depths][::-1]
legend_lines = legend_lines[::-1]  # Reverse the handles to match label order

fig.legend(
    legend_lines,
    legend_labels,
    loc="center left",
    bbox_to_anchor=(1.005, 0.5),
    fontsize=12,
    title="Read depth",
    title_fontsize=12,
    frameon=False
)
# Save in PDF and PNG formats
output_base = "simulation_results/simulated_read_lengths_all"
plt.savefig(f"{output_base}.pdf", bbox_inches='tight')
plt.savefig(f"{output_base}.png", dpi=600, bbox_inches='tight')
plt.close(fig)
import matplotlib.pyplot as plt
import numpy as np
import json
import math
import matplotlib.cm as cm
from matplotlib import gridspec
import seaborn as sns

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 30})

# Load the data
with open("/hps/nobackup/iqbal/dander/thesis_figures/sim_evaluation/simulation_evaluation/simulation_results/copy_number_tuples.json") as i:
    data = json.load(i)

# Filter out scenario "0"
data = {k: v for k, v in data.items() if k != "1"}

# Set up the figure with subplots for both scatter and KDE
num_scenarios = len(data)
fig = plt.figure(figsize=(13.5, 45))  # slightly less wide
gs = gridspec.GridSpec(num_scenarios, 2, width_ratios=[2, 1])

# Loop over scenarios
for idx, (scenario, depths) in enumerate(data.items()):
    scatter_ax = fig.add_subplot(gs[idx, 0])
    kde_ax = fig.add_subplot(gs[idx, 1])
    colors = sns.color_palette("colorblind", len(depths))
    max_pred = 0
    plotted = False

    for color, (depth_str, pairs) in zip(colors, sorted(depths.items(), key=lambda x: int(x[0]))):
        if not pairs:
            continue
        
        depth = int(depth_str)
        true_vals = np.array([p[0] for p in pairs])
        pred_vals = np.array([p[1] for p in pairs])
        max_pred = max(max_pred, max(pred_vals, default=0))

        # Scatter plot
        scatter_ax.scatter(true_vals, pred_vals, label=f'Depth {depth}', alpha=1, s=50, color=color)
        plotted = True

        # Seaborn KDE plot
        if len(pred_vals) > 1:
            sns.kdeplot(
                y=pred_vals,
                ax=kde_ax,
                fill=True,
                alpha=0.9,
                linewidth=0,
                color=color,
                label=f'Depth {depth}',
                common_norm=False,
                bw_adjust=0.3  # Optional: adjust bandwidth
            )

    # Formatting scatter plot
    scatter_ax.set_title(f"Scenario {int(scenario) - 1}", fontsize=30)
    if idx == 4:
        scatter_ax.set_xlabel("True cellular copy number")
    else:
        scatter_ax.set_xlabel("")
    scatter_ax.set_ylabel("Amira cellular copy number estimate")
    scatter_ax.set_xlim([0, math.ceil(max_pred) + 1])
    scatter_ax.set_ylim([0, math.ceil(max_pred) + 1])
    scatter_ax.spines['top'].set_visible(False)
    scatter_ax.spines['right'].set_visible(False)
    #if plotted:
     #   scatter_ax.legend()
    scatter_ax.grid(True)

    # Format KDE plot: remove all axis elements
    kde_ax.set_ylim([0, math.ceil(max_pred) + 1])
    kde_ax.set_xticks([])
    kde_ax.set_yticks([])
    kde_ax.set_title("")
    kde_ax.set_xlabel("")
    for spine in kde_ax.spines.values():
        spine.set_visible(False)
    kde_ax.yaxis.grid(True)
    kde_ax.xaxis.grid(False)
    kde_ax.legend().set_visible(False)

from matplotlib.lines import Line2D

# Collect legend handles/labels from the last set of colors used
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label=f'Depth {int(depth)}',
           markerfacecolor=color, markersize=8)
    for depth, color in zip(sorted([int(d) for d in depths.keys()]), cm.viridis(np.linspace(0, 1, len(depths))))
]

plt.tight_layout()  # Leave space on the right
plt.savefig("/hps/nobackup/iqbal/dander/thesis_figures/sim_evaluation/simulation_evaluation/simulation_results/copy_numbers_by_scenario.png", dpi=600)
plt.savefig("/hps/nobackup/iqbal/dander/thesis_figures/sim_evaluation/simulation_evaluation/simulation_results/copy_numbers_by_scenario.pdf")
plt.close()

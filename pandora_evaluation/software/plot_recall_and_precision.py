import os
import matplotlib.pyplot as plt
import re
import pandas as pd
import glob
from tqdm import tqdm
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


# Directory containing the files
input_directory = glob.glob("pandora_assessment_with_AMR_filtered/*/summary_stats.k*.w*.txt")
# Regex to extract values of w, k, per gene recall, and per gene precision
pattern = re.compile(r"Per gene recall: ([0-9.]+).*Per gene precision: ([0-9.]+)", re.DOTALL)

# Initialize a list to store the extracted data
data = []

# Loop through each file in the directory
for file_path in tqdm(input_directory):
    with open(file_path, "r") as file:
        content = file.read()
        match = pattern.search(content)
        if match:
            w = int(os.path.basename(file_path).split(".w")[1].split(".txt")[0])
            k = int(os.path.basename(file_path).split(".k")[1].split(".w")[0])
            c = float(os.path.basename(os.path.dirname(file_path)).split("_")[0])
            l = float(os.path.basename(os.path.dirname(file_path)).split("_")[1])
            recall = float(match.group(1))
            precision = float(match.group(2))
            data.append({"w": w, "k": k, "recall": recall, "precision": precision, "identity": c, "length": l})

data = sorted(data, key=lambda x: (x["w"], x["k"]))
# Convert the data into a pandas DataFrame
df = pd.DataFrame(data)
# Set global font properties
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 18})

# First figure: Recall and precision vs. length for identity groups
fig1, ax1 = plt.subplots(1, 2, figsize=(12, 8))  # Square subplots

# Recall (left subplot)
max_val = 0
max_c = 0
max_l = 0
min_val = 10000000
min_c = 0
min_l = 0
for i_value, group in df.groupby("identity"):
    filtered_group = group[(group["w"] == 14) & (group["k"] == 15)]
    if not filtered_group.empty:
        ax1[0].plot(filtered_group["length"], filtered_group["recall"], label=f"Cluster identity={i_value:.2f}", marker="o", markersize=8)
        # Find max recall in the filtered group
        max_recall = filtered_group["recall"].max()
        max_recall_row = filtered_group[filtered_group["recall"] == max_recall].iloc[0]
        if max_recall > max_val:
            max_val = max_recall
            max_c = i_value
            max_l = max_recall_row["length"]

        # Find min recall in the filtered group
        min_recall = filtered_group["recall"].min()
        min_recall_row = filtered_group[filtered_group["recall"] == min_recall].iloc[0]
        if min_recall < min_val:
            min_val = min_recall
            min_c = i_value
            min_l = min_recall_row["length"]
print(f"\nMax recall at c={max_c} l={max_l}: {max_val}, min recall at c={min_c} l={min_l}: {min_val} \n")
ax1[0].set_xlabel("Cluster length threshold", fontsize=18)
ax1[0].set_ylabel("Recall", fontsize=18)
ax1[0].grid(True)
ax1[0].tick_params(axis='both', labelsize=18)
ax1[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
ax1[0].set_box_aspect(1)
ax1[0].yaxis.set_major_locator(MultipleLocator(0.05))  # Set y-tick every 0.05
ax1[0].set_ylim([0.8, 1])
ax1[0].set_xlim([-0.05, 1])
ax1[1].xaxis.set_major_locator(MultipleLocator(0.2))

# Precision (right subplot)
max_val = 0
max_c = 0
max_l = 0
min_val = 10000000
min_c = 0
min_l = 0
for i_value, group in df.groupby("identity"):
    filtered_group = group[(group["w"] == 14) & (group["k"] == 15)]
    if not filtered_group.empty:
        ax1[1].plot(filtered_group["length"], filtered_group["precision"], marker="o", markersize=8)
        # Find max recall in the filtered group
        max_precision = filtered_group["precision"].max()
        max_precision_row = filtered_group[filtered_group["precision"] == max_precision].iloc[0]
        if max_precision > max_val:
            max_val = max_precision
            max_c = i_value
            max_l = max_precision_row["length"]

        # Find min recall in the filtered group
        min_precision = filtered_group["precision"].min()
        min_precision_row = filtered_group[filtered_group["precision"] == min_precision].iloc[0]
        if min_precision < min_val:
            min_val = min_precision
            min_c = i_value
            min_l = min_precision_row["length"]
print(f"\nMax precision at c={max_c} l={max_l}: {max_val}, min precision at c={min_c} l={min_l}: {min_val} \n")
ax1[1].set_xlabel("Cluster length threshold", fontsize=18)
ax1[1].set_ylabel("Precision", fontsize=18)
ax1[1].grid(True)
ax1[1].tick_params(axis='both', labelsize=18)
ax1[1].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
ax1[1].set_box_aspect(1)
ax1[1].yaxis.set_major_locator(MultipleLocator(0.05))  # Set y-tick every 0.05
ax1[1].set_ylim([0.8, 1])
ax1[1].set_xlim([-0.05, 1])
ax1[1].xaxis.set_major_locator(MultipleLocator(0.2))

# Add a single legend for the entire figure
fig1.legend(loc="center left", bbox_to_anchor=(1.005, 0.5), fontsize=18, title="Cluster Identity", title_fontsize=18)
fig1.tight_layout()
fig1.savefig("result_plots_identity_vs_length.png", dpi=600, bbox_inches="tight")
fig1.savefig("result_plots_identity_vs_length.pdf", bbox_inches="tight")
plt.close()

# Second figure: Recall and precision vs. w for k groups
fig2, ax2 = plt.subplots(1, 2, figsize=(12, 8))  # Square subplots

# Recall (left subplot)
max_val = 0
max_c = 0
max_l = 0
min_val = 10000000
min_c = 0
min_l = 0
for k_value, group in df.groupby("k"):
    filtered_group = group[(group["identity"] == 0.80) & (group["length"] == 0.0)]
    if not filtered_group.empty:
        ax2[0].plot(filtered_group["w"], filtered_group["recall"], label=f"$k$={k_value}", marker="o", markersize=8)
        # Find max recall in the filtered group
        max_recall = filtered_group["recall"].max()
        max_recall_row = filtered_group[filtered_group["recall"] == max_recall].iloc[0]
        if max_recall > max_val:
            max_val = max_recall
            max_k = k_value
            max_w = max_recall_row["w"]

        # Find min recall in the filtered group
        min_recall = filtered_group["recall"].min()
        min_recall_row = filtered_group[filtered_group["recall"] == min_recall].iloc[0]
        if min_recall < min_val:
            min_val = min_recall
            min_k = k_value
            min_w = min_recall_row["w"]
print(f"\nMax recall at k={max_k} w={max_w}: {max_val}, min recall at k={min_k} w={min_w}: {min_val} \n")
ax2[0].set_xlabel("$w$", fontsize=18)
ax2[0].set_ylabel("Recall", fontsize=18)
ax2[0].grid(True)
ax2[0].tick_params(axis='both', labelsize=18)
ax2[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
ax2[0].set_box_aspect(1)
ax2[0].set_ylim([0.6, 1])

# Precision (right subplot)
max_val = 0
max_c = 0
max_l = 0
min_val = 10000000
min_c = 0
min_l = 0
for k_value, group in df.groupby("k"):
    filtered_group = group[(group["identity"] == 0.80) & (group["length"] == 0.0)]
    if not filtered_group.empty:
        ax2[1].plot(filtered_group["w"], filtered_group["precision"], marker="o", markersize=8)
        # Find max recall in the filtered group
        max_precision = filtered_group["precision"].max()
        max_precision_row = filtered_group[filtered_group["precision"] == max_precision].iloc[0]
        if max_precision > max_val:
            max_val = max_precision
            max_k = k_value
            max_w = max_precision_row["w"]

        # Find min recall in the filtered group
        min_precision = filtered_group["precision"].min()
        min_precision_row = filtered_group[filtered_group["precision"] == min_precision].iloc[0]
        if min_precision < min_val:
            min_val = min_precision
            min_k = k_value
            min_w = min_precision_row["w"]
print(f"\nMax precision at k={max_k} w={max_w}: {max_val}, min precision at k={min_k} w={min_w}: {min_val} \n")
ax2[1].set_xlabel("$w$", fontsize=18)
ax2[1].set_ylabel("Precision", fontsize=18)
ax2[1].grid(True)
ax2[1].tick_params(axis='both', labelsize=18)
ax2[1].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
ax2[1].set_box_aspect(1)
ax2[1].set_ylim([0.6, 1])

# Add a single legend for the entire figure
fig2.legend(loc="center left", bbox_to_anchor=(1.005, 0.5), fontsize=18, title="$k$-mer Size", title_fontsize=18)
fig2.tight_layout()
fig2.savefig("result_plots_k_vs_w.png", dpi=600, bbox_inches="tight")
fig2.savefig("result_plots_k_vs_w.pdf", bbox_inches="tight")
plt.close()

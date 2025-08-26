import os
import pandas as pd
import matplotlib.pyplot as plt
import json

amira = "amira_output"
flye = "AMRFinderPlus_flye"
metaflye = "AMRFinderPlus_metaflye"
nanomdbg = "AMRFinderPlus_nanomdbg"
truth = "truth_jsons/RW_data/all_reads.json"

def apply_rules(gene):
    if "blaCTX-M" in gene:
        gene = "blaCTX-M"
    if "blaNDM" in gene:
        gene = "blaNDM"
    if "blaOXA" in gene:
        gene = "blaOXA"
    if "aac(6')-Ib" in gene:
        gene = "aac(6')-Ib"
    if "blaEC" in gene:
        gene = "blaEC"
    if "oqx" in gene:
        gene = "oqx"
    if "blaTEM" in gene:
        gene = "blaTEM"
    if "fosA" in gene:
        gene = "fosA"
    if "blaCMY" in gene:
        gene = "blaCMY"
    if "aadA" in gene:
        gene = "aadA"
    if "arr-" in gene:
        gene = "arr-"
    if "dfrA" in gene:
        gene = "dfrA"
    if "rmtB" in gene:
        gene = "rmtB"
    if "aac(3)-II" in gene and "aac(3)-III" not in gene:
        gene = "aac(3)-II"
    if "aac(3)-III" in gene:
        gene = "aac(3)-III"
    if "blaSHV" in gene:
        gene = "blaSHV"
    if "qnrS" in gene:
        gene = "qnrS"
    if "aac(3)-I" in gene and "aac(3)-II" not in gene and "aac(3)-III" not in gene:
        gene = "aac(3)-I"
    if "blaKPC" in gene:
        gene = "blaKPC"
    if "mcr-5" in gene:
        gene = "mcr-5"
    if "qnrB" in gene:
        gene = "qnrB"
    if "cmlA" in gene:
        gene = "cmlA"
    if "aph(3'')-I" in gene and "aph(3'')-II" not in gene and "aph(3'')-III" not in gene:
        gene = "aph(3)-I"
    if "aph(3'')-II" in gene and "aph(3'')-III" not in gene:
        gene = "aph(3)-II"
    if "aph(3'')-III" in gene:
        gene = "aph(3)-III"
    if "aph(3')-I" in gene and "aph(3')-II" not in gene and "aph(3')-III" not in gene:
        gene = "aph(3)-I"
    if "aph(3')-II" in gene and "aph(3')-III" not in gene:
        gene = "aph(3)-II"
    if "aph(3')-III" in gene:
        gene = "aph(3)-III"
    if "aph(4)-I" in gene and "aph(4)-II" not in gene and "aph(4)-III" not in gene:
        gene = "aph(4)-I"
    if "aph(4)-II" in gene and "aph(4)-III" not in gene:
        gene = "aph(4)-II"
    if "aph(4)-III" in gene:
        gene = "aph(4)-III"
    if "aph(6)-I" in gene and "aph(6)-II" not in gene and "aph(6)-III" not in gene:
        gene = "aph(6)-I"
    if "aph(6)-II" in gene and "aph(6)-III" not in gene:
        gene = "aph(6)-II"
    if "aph(6)-III" in gene:
        gene = "aph(6)-III"
    if "qacE" in gene:
        gene = "qacE"
    if "blaLAP" in gene:
        gene = "blaLAP"
    if "aac(6')-I" in gene and "aac(6')-II" not in gene and "aac(6')-III" not in gene:
        gene = "aac(6')-I"
    if "aac(6')-II" in gene and "aac(6')-III" not in gene:
        gene = "aac(6')-II"
    if "aac(6')-III" in gene:
        gene = "aac(6')-III"
    if "blaDHA" in gene:
        gene = "blaDHA"
    if "qepA" in gene:
        gene = "qepA"
    if "blaIMI" in gene:
        gene = "blaIMI"
    return gene
# Load truth data
truth_counts = {}
with open(truth) as i:
    truth = json.load(i)
for g in truth:
    truth_counts[apply_rules(g)] = truth[g]

depths = [10, 20, 40, 60, 80, 100, 120, 140]
amira_recalls, mdbg_recalls, mflye_recalls, flye_recalls = [], [], [], []
amira_precisions, mdbg_precisions, mflye_precisions, flye_precisions = [], [], [], []

for depth in depths:
    def load_counts(path, gene_col, coverage_col=None):
        counts = {}
        if os.path.exists(path):
            data = pd.read_csv(path, sep="\t")
            for _, row in data.iterrows():
                if coverage_col and (row[coverage_col] < 90 or row["Element subtype"] == "POINT"):
                    continue
                gene = apply_rules(row[gene_col])
                counts[gene] = counts.get(gene, 0) + 1
        return counts

    flye_counts = load_counts(os.path.join(flye, f"RW_data.{depth}x_reads"), "Gene symbol", "% Coverage of reference sequence")
    mflye_counts = load_counts(os.path.join(metaflye, f"RW_data.{depth}x_reads"), "Gene symbol", "% Coverage of reference sequence")
    mdbg_counts = load_counts(os.path.join(nanomdbg, f"RW_data.{depth}x_reads"), "Gene symbol", "% Coverage of reference sequence")
    amira_counts = load_counts(os.path.join(amira, f"RW_data.{depth}x_reads", "amira_results.tsv"), "Determinant name")

    def compute_metrics(pred_counts):
        tp = sum(min(pred_counts.get(g, 0), truth_counts[g]) for g in truth_counts)
        fp = sum(
            (pred_counts.get(g, 0) - truth_counts[g])
            for g in truth_counts
            if pred_counts.get(g, 0) > truth_counts[g]
        )
        recall = tp / sum(truth_counts.values()) * 100
        precision = tp / (tp + fp) * 100 if (tp + fp) > 0 else 0
        return recall, precision

    ar, ap = compute_metrics(amira_counts)
    mr, mp = compute_metrics(mdbg_counts)
    mfr, mfp = compute_metrics(mflye_counts)
    fr, fp = compute_metrics(flye_counts)

    amira_recalls.append(ar)
    amira_precisions.append(ap)
    mdbg_recalls.append(mr)
    mdbg_precisions.append(mp)
    mflye_recalls.append(mfr)
    mflye_precisions.append(mfp)
    if depth < 80:
        flye_recalls.append(fr)
        flye_precisions.append(fp)

for i in range(len(depths)):
    try:
        print(f"Depth: {depths[i]}x")
        print(f"Amira recall: {amira_recalls[i]}, Flye recall: {flye_recalls[i]}, metaFlye recall: {mflye_recalls[i]}, nanoMDBG recall: {mdbg_recalls[i]}")
        print(f"Amira precision: {amira_precisions[i]}, Flye precision: {flye_precisions[i]}, metaFlye precision: {mflye_precisions[i]}, nanoMDBG precision: {mdbg_precisions[i]}\n")
    except:
        print(f"Amira recall: {amira_recalls[i]}, metaFlye recall: {mflye_recalls[i]}, nanoMDBG recall: {mdbg_recalls[i]}")
        print(f"Amira precision: {amira_precisions[i]}, metaFlye precision: {mflye_precisions[i]}, nanoMDBG precision: {mdbg_precisions[i]}\n")

# Plotting
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 16})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 10), sharex=True)

colors = ['#0173b2', '#de8f05', '#029e73', '#d55e00']
# Recall plot
l1, = ax1.plot(depths, amira_recalls, label='Amira', marker='o', linewidth=1.5, color=colors[0])
l2, = ax1.plot(depths[:4], flye_recalls, label='Flye', marker='X', linewidth=1.5, color=colors[3])
l3, = ax1.plot(depths, mflye_recalls, label='metaFlye', marker='^', linewidth=1.5, color=colors[1])
l4, = ax1.plot(depths, mdbg_recalls, label='nanoMDBG', marker='s', linewidth=1.5, color=colors[2])
ax1.set_xlabel('')
ax1.set_ylabel('Recall (%)')
ax1.set_xlim([0, 145])
ax1.set_ylim([0, 105])
ax1.set_xticks([0, 20, 40, 60, 80, 100, 120, 140])
ax1.set_yticks([0, 20, 40, 60, 80, 100])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.tick_params(axis='both', which='major', length=6, width=1)

# Precision plot
ax2.plot(depths, amira_precisions, label='Amira', marker='o', linewidth=1.5)
ax2.plot(depths[:4], flye_precisions, label='Flye', marker='X', linewidth=1.5)
ax2.plot(depths, mflye_precisions, label='metaFlye', marker='^', linewidth=1.5)
ax2.plot(depths, mdbg_precisions, label='nanoMDBG', marker='s', linewidth=1.5)
ax2.set_xlabel('Sequencing depth (×)')
ax2.set_ylabel('Precision (%)')
ax2.set_xlim([0, 145])
ax2.set_ylim([0, 105])
ax2.set_xticks([0, 20, 40, 60, 80, 100, 120, 140])
ax2.set_yticks([0, 20, 40, 60, 80, 100])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.tick_params(axis='both', which='major', length=6, width=1)

# Single legend to the right
fig.legend([l1, l2, l3, l4], ['Amira', 'Flye', 'metaFlye', 'nanoMDBG'],
           loc='center left', bbox_to_anchor=(0.9, 0.5), frameon=False)
#fig.legend([l1, l3, l4], ['Amira', 'metaFlye', 'nanoMDBG'],
#           loc='center left', bbox_to_anchor=(0.9, 0.5), frameon=False)
plt.savefig("summary_plot_precision_recall.png", dpi=900, bbox_inches='tight')
plt.savefig("summary_plot_precision_recall.pdf", bbox_inches='tight')

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 14,
    'axes.linewidth': 1,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.frameon': False,
})

def process_AMRFP_results(file_list, reference_genes):
    merged_results = {}
    for f in file_list:
        sample = os.path.basename(f).replace(".tsv", "")
        merged_results[sample] = {}
        results = pd.read_csv(f, sep="\t")
        for index, row in results.iterrows():
            gene = row["Gene symbol"]
            if row["Element subtype"] == "POINT":
                continue
            if "PARTIAL" in row["Method"]:
                continue
            gene = apply_rules(gene)
            if gene not in reference_genes:
                continue
            if gene not in merged_results[sample]:
                merged_results[sample][gene] = 0
            merged_results[sample][gene] += 1
    return merged_results

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
    if "oqx" in gene or "Oqx" in gene:
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
    if "ant(3'')-I" in gene and "ant(3'')-II" not in gene and "ant(3'')-III" not in gene:
        gene = "ant(3'')-I"
    if "ant(3'')-II" in gene and "ant(3'')-III" not in gene:
        gene = "ant(3'')-II"
    if "ant(3'')-III" in gene:
        gene = "ant(3'')-III"
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
    if "ant(6)-I" in gene and "ant(6)-II" not in gene and "ant(6)-III" not in gene:
        gene = "ant(6)-I"
    if "ant(6)-II" in gene and "ant(6)-III" not in gene:
        gene = "ant(6)-II"
    if "ant(6)-III" in gene:
        gene = "ant(6)-III"
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
    if "blaPDC" in gene:
        gene = "blaPDC"
    if "blaADC" in gene:
        gene = "blaADC"
    if "blaPC" in gene:
        gene = "blaPC"
    if "blaR" in gene:
        gene = "blaR"
    if "catA" in gene:
        gene = "catA"
    if "catB" in gene:
        gene = "catB"
    if "mcr-1" in gene:
        gene = "mcr-1"
    return gene

def load_amira_results(amira_dirs):
    amira_results = {}
    for s in amira_dirs:
        if os.path.exists(os.path.join(s, "amira_results.tsv")):
            sample = os.path.basename(s.replace("amira_output.", ""))
            amira_results[sample] = {}
            amira_table = pd.read_csv(os.path.join(s, "amira_results.tsv"), sep="\t")
            for index, row in amira_table.iterrows():
                reference_gene = apply_rules(row["Determinant name"])
                if reference_gene not in amira_results[sample]:
                    amira_results[sample][reference_gene] = 0
                amira_results[sample][reference_gene] += 1
    return amira_results

def calculate_concordance(nanomdbg_results, metaflye_results, amira_results, reference_genes):
    stats = {
        "nanomdbg": {"presence_agree": 0, "copy_agree": 0},
        "metaflye": {"presence_agree": 0, "copy_agree": 0},
        "amira": {"presence_agree": 0, "copy_agree": 0},
    }

    all_samples = sorted(set(nanomdbg_results) & set(metaflye_results) & set(amira_results))
    nano_genes = set()
    meta_genes = set()
    amira_genes = set()
    for s in all_samples:
        nano_genes.update(list(nanomdbg_results[sample].keys()))
        meta_genes.update(list(metaflye_results[sample].keys()))
        amira_genes.update(list(amira_results[sample].keys()))

    for g in nano_genes:
        if all(g in nanomdbg_results[s] for s in all_samples):
            stats["nanomdbg"]["presence_agree"] += 1
        if len(set([nanomdbg_results[s].get(g, 0) for s in all_samples])) == 1:
            stats["nanomdbg"]["copy_agree"] += 1
    for g in meta_genes:
        if all(g in metaflye_results[s] for s in all_samples):
            stats["metaflye"]["presence_agree"] += 1
        if len(set([metaflye_results[s].get(g, 0) for s in all_samples])) == 1:
            stats["metaflye"]["copy_agree"] += 1
    for g in amira_genes:
        if all(g in amira_results[s] for s in all_samples):
            stats["amira"]["presence_agree"] += 1
        if len(set([amira_results[s].get(g, 0) for s in all_samples])) == 1:
            stats["amira"]["copy_agree"] += 1

    print("\n### Method-wise Concordance Across All Samples ###")
    for method in stats:
        if method == "nanomdbg":
            total_pairs = len(nano_genes)
        if method == "metaflye":
            total_pairs = len(meta_genes)
        if method == "amira":
            total_pairs = len(amira_genes)
        pres_prop = stats[method]["presence_agree"] / total_pairs if total_pairs > 0 else 0
        copy_prop = stats[method]["copy_agree"] / total_pairs if total_pairs > 0 else 0
        print(f"Method: {method}")
        print(f"  Total genes: {total_pairs}")
        print(f"  Presence/Absence agreement: {stats[method]['presence_agree']} ({pres_prop:.2%})")
        print(f"  Copy number agreement:      {stats[method]['copy_agree']} ({copy_prop:.2%})\n")



# File paths
allele_file = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/metagenome_panRG_thesis/AMR_alleles_unified.fa"
metaflye_dir = "/hps/nobackup/iqbal/dander/amira_metagenome/Purushothaman_data/metaflye_amr"
nanomdbg_dir = "/hps/nobackup/iqbal/dander/amira_metagenome/Purushothaman_data/nanomdbg_amr"

# Load reference genes
with open(allele_file) as i:
    allele_rows = i.read().split(">")[1:]
reference_genes = set()
for r in allele_rows:
    if r != "":
        amira_allele, reference_allele = r.split("\n")[0].split(";")
        reference_genes.add(apply_rules(reference_allele.split(".NG")[0]))

# Load tool results
metaflye_results = process_AMRFP_results(glob.glob(os.path.join(metaflye_dir, "*1.tsv")), reference_genes)
nanomdbg_results = process_AMRFP_results(glob.glob(os.path.join(nanomdbg_dir, "*1.tsv")), reference_genes)
amira_results = load_amira_results(glob.glob("/hps/nobackup/iqbal/dander/amira_metagenome/Purushothaman_data/amira_output*1"))

# Combine into multi-panel figure
all_samples = sorted(set(metaflye_results.keys()) | set(nanomdbg_results.keys()) | set(amira_results.keys()))
num_samples = len(all_samples)

fig, axes = plt.subplots(num_samples, 1, figsize=(14, 5 * num_samples))

if num_samples == 1:
    axes = [axes]  # Make iterable if only one sample

legend_handles, legend_labels = None, None

color_palette = ListedColormap(['#0173b2', '#de8f05', '#029e73'])  # blue, orange, green

for idx, sample in enumerate(all_samples):
    data = []
    for gene in sorted(reference_genes):
        amira_count = amira_results.get(sample, {}).get(gene, 0)
        metaflye_count = metaflye_results.get(sample, {}).get(gene, 0)
        nanomdbg_count = nanomdbg_results.get(sample, {}).get(gene, 0)
        if nanomdbg_count + metaflye_count + amira_count > 0:
            data.append([gene, amira_count, metaflye_count, nanomdbg_count])

    if not data:
        continue

    df = pd.DataFrame(data, columns=["Gene", "Amira", "metaFlye", "nanoMDBG"])
    df.set_index("Gene", inplace=True)

    ax = axes[idx]
    df.plot(kind="bar", stacked=True, ax=ax, width=0.8, color=color_palette.colors)
    ax.set_title(f"Sample: {sample}")
    ax.set_xticklabels([f"$\\it{{{label}}}$" for label in df.index], rotation=90)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(False)
    if idx == 3:
        ax.set_xlabel("Gene")
    else:
        ax.set_xlabel("")
    
    if legend_handles is None and legend_labels is None:
        legend_handles, legend_labels = ax.get_legend_handles_labels()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    ax.get_legend().remove()

# Add single legend on the right center
fig.legend(legend_handles, legend_labels, loc="center right", bbox_to_anchor=(1, 0.5), frameon=False)
fig.supylabel("Genomic copy number count")

plt.tight_layout(rect=[0, 0, 0.87, 1])
plt.savefig("amr_comparison_multi_panel_stacked.png")
plt.savefig("amr_comparison_multi_panel_stacked.pdf")
plt.close()

calculate_concordance(nanomdbg_results, metaflye_results, amira_results, reference_genes)
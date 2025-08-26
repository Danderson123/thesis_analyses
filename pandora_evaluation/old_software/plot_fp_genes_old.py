import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import numpy as np
from matplotlib import gridspec

# Parameters
gene = "group_160"
kmer_size = 15  # k-mer size
flank = 0  # flanking region around hit positions
max_read_len = 2000  # maximum length to visualize

# File paths
base_dir = "/hps/nobackup/iqbal/dander/thesis_figures/pandora_systematic_evaluation/pandora_assessment_with_AMR/0.8_0/ad_hoc"
minimatches_file = f"{base_dir}/{gene}_pandora/Escherichia_coli_MINF_1A.{gene}.fastq.minimatches"
reads_file = f"{base_dir}/Escherichia_coli_MINF_1A.{gene}.fastq.gz"
consensus_file = f"{base_dir}/{gene}_pandora/pandora.consensus.fq.gz"
outfile_base = f"{base_dir}/{gene}_hits"
fasta_file = f"/hps/nobackup/iqbal/dander/thesis_figures/pandora_systematic_evaluation/pandora_assessment_with_AMR/0.8_0/panaroo_output/aligned_gene_sequences/{gene}.fasta"


# Read minimatches and collect hits
mm_df = pd.read_csv(minimatches_file, sep="\t")
hits = {}
positions = {}
all_kmers = set()

def rc(seq):
    replacement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    seq_list = reversed(list(seq))
    return "".join([replacement[b] for b in seq_list])

for index, row in tqdm(mm_df.iterrows(), total=len(mm_df)):
    if row["prg"] == gene and len(hits) < 50:
        read = row["read"]
        if read not in hits:
            hits[read] = set()
            positions[read] = {"start": float('inf'), "end": 0}
        hits[read].add(row["kmer"])
        all_kmers.add(row["kmer"])
        all_kmers.add(rc(row["kmer"]))
        positions[read]["start"] = min(positions[read]["start"], row["read_start"])
        positions[read]["end"] = max(positions[read]["end"], row["read_end"])
    elif len(hits) >= 50:
        break

print(f"Number of selected reads: {len(hits)}")

# Load selected read sequences
fastq_dict = {}
with pysam.FastxFile(reads_file) as fh:
    for entry in fh:
        if entry.name in hits:
            fastq_dict[entry.name] = {"sequence": entry.sequence, "quality": entry.quality}

fasta_content = pysam.FastaFile(fasta_file)
fasta_dict = {}
hits_ref = {}

for ref in fasta_content.references:
    seq = fasta_content[ref].replace("-", "")
    fasta_dict[ref] = {"sequence": seq}
    hits_ref[ref] = set()

# Prepare plot data
read_order = []
matched_points = {}  # kmer -> list of (x, y)
unmatched_x, unmatched_y = [], []

for read_name in reversed(list(fastq_dict.keys())):
    seq = fastq_dict[read_name]["sequence"]
    start = max(0, positions[read_name]["start"] - flank)
    end = min(len(seq), positions[read_name]["end"] + flank)
    region = seq[start:end]
    if len(region) > max_read_len:
        continue
    kmers = [region[j:j+kmer_size] for j in range(len(region) - kmer_size + 1)]
    row_matched = False
    for j, kmer in enumerate(kmers):
        if kmer in hits[read_name]:
            matched_points.setdefault(kmer, []).append((j, len(read_order)))
            row_matched = True
        else:
            unmatched_x.append(j)
            unmatched_y.append(len(read_order))
    if row_matched or any(kmer not in hits[read_name] for kmer in kmers):
        read_order.append(read_name)

fasta_content = pysam.FastaFile(fasta_file)
fasta_dict = {}
ref_order = []
matched_points_ref = {}  # kmer -> list of (x, y)
unmatched_x_ref, unmatched_y_ref = [], []
for ref in fasta_content.references:
   # if len(ref_order)< 10:
    sequence = fasta_content[ref].replace("-", "").upper()
    kmers = [sequence[j:j+kmer_size] for j in range(len(sequence) - kmer_size + 1)]
    if any(k in all_kmers for k in kmers):
        for j, kmer in enumerate(kmers):
            if kmer in all_kmers:
                matched_points_ref.setdefault(kmer, []).append((j, len(ref_order)))
            elif rc(kmer) in all_kmers:
                matched_points_ref.setdefault(rc(kmer), []).append((j, len(ref_order)))
            else:
                unmatched_x_ref.append(j)
                unmatched_y_ref.append(len(ref_order))
        ref_order.append(ref)

# Plotting
fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
ax0 = plt.subplot(gs[0])  # reference
ax1 = plt.subplot(gs[1])  # reads

# Reference panel
ax0.scatter(unmatched_x_ref, unmatched_y_ref, c="lightgrey", label="Unmatched", s=12, alpha=0.3)
colormap = plt.cm.get_cmap("tab10", len(matched_points_ref))
for idx, (kmer, coords) in enumerate(matched_points_ref.items()):
    x, y = zip(*coords)
    ax0.scatter(x, y, s=12, alpha=1, color=colormap(idx))
ax0.set_title("$k$-mer hits to references")
ax0.set_ylabel("reference")
ax0.set_yticks([])
ax0.set_yticklabels([])
ax0.set_xlim([0, 300])
ax0.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)


ax1.scatter(unmatched_x, unmatched_y, c="lightgrey", s=12, label="Unmatched", alpha=0.3)
# Assign a unique color per k-mer
colormap = plt.cm.get_cmap("tab10", len(matched_points))
for idx, (kmer, coords) in enumerate(matched_points.items()):
    x, y = zip(*coords)
    ax1.scatter(x, y, s=12, alpha=1, label=kmer, color=colormap(idx))

ax1.set_xlabel("$k$-mer index")
ax1.set_ylabel("read")
ax1.set_title("$k$-mer hits for false positive calls in the reads")

ax1.set_yticks([])
ax1.set_yticklabels([])
ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
ax1.set_xlim([0, 300])
plt.tight_layout()
# Save outputs
plt.savefig(f"{outfile_base}.png", dpi=600)
plt.savefig(f"{outfile_base}.pdf")
print(f"Plots saved to {outfile_base}.png/.pdf")

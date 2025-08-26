import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import numpy as np
from matplotlib import gridspec

# Parameters
gene = "group_227"
kmer_size = 15  # k-mer size
flank = 0  # flanking region around hit positions
max_read_len = 2000  # maximum length to visualize

# File paths
base_dir = "/hps/nobackup/iqbal/dander/thesis_figures/pandora_systematic_evaluation/pandora_assessment_with_AMR/0.8_0/ad_hoc"
minimatches_file = f"{base_dir}/{gene}_pandora/Escherichia_coli_MINF_1A.{gene}.fastq.minimatches"
minimatches_file_tp = f"{base_dir}/{gene}_tp_pandora/{gene}_tp_reads.fastq.minimatches"
reads_file = f"{base_dir}/Escherichia_coli_MINF_1A.{gene}.fastq.gz"
tp_reads_file = f"{base_dir}/{gene}_tp_reads.fastq.gz"
consensus_file = f"{base_dir}/{gene}_pandora/pandora.consensus.fq.gz"
outfile_base = f"{base_dir}/{gene}_hits"
fasta_file = f"/hps/nobackup/iqbal/dander/thesis_figures/pandora_systematic_evaluation/pandora_assessment_with_AMR/0.8_0/panaroo_output/aligned_gene_sequences/{gene}.fasta"

# Reverse complement function
def rc(seq):
    replacement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join([replacement[b] for b in reversed(seq)])

def load_fasta_and_get_kmers(fasta_file, kmer_size):
    # Load reference sequences and extract all kmers
    fasta_content = pysam.FastaFile(fasta_file)
    fasta_dict = {}
    all_kmers = set()
    for ref in fasta_content.references:
        sequence = fasta_content[ref].replace("-", "").upper()
        fasta_dict[ref] = sequence
        kmers = [sequence[j:j+kmer_size] for j in range(len(sequence) - kmer_size + 1)]
        for kmer in kmers:
            all_kmers.add(kmer)
            all_kmers.add(rc(kmer))
    return fasta_dict, all_kmers, fasta_dict

def load_fastq_and_get_kmers(minimatches_file, reads_file, hit_kmers = set()):
    # Load minimatches
    mm_df = pd.read_csv(minimatches_file, sep="\t")
    # Collect hits from minimatches
    hits = {}
    positions = {}
    for _, row in tqdm(mm_df.iterrows(), total=len(mm_df)):
        if row["prg"] == gene and len(hits) < 50:
            read = row["read"]
            if read not in hits:
                hits[read] = set()
                positions[read] = {"start": float('inf'), "end": 0}
            if row["kmer"] in all_kmers:
                positions[read]["start"] = min(positions[read]["start"], row["read_start"])
                positions[read]["end"] = max(positions[read]["end"], row["read_end"])
                hits[read].add(row["kmer"])
                hit_kmers.add(row["kmer"])
        elif len(hits) >= 50:
            break

    # Load selected read sequences
    fastq_dict = {}
    with pysam.FastxFile(reads_file) as fh:
        for entry in fh:
            if entry.name in hits:
                fastq_dict[entry.name] = {"sequence": entry.sequence, "quality": entry.quality}

    # Read hits: collect matched/unmatched coords + matched_points dict
    unmatched_x, unmatched_y = [], []
    matched_x, matched_y = [], []
    matched_points = {}  # kmer -> list of (x, y)
    read_order = list(hits.keys())

    count = 0
    for i, read_name in enumerate(read_order):
        seq = fastq_dict[read_name]["sequence"]
        start = max(0, positions[read_name]["start"] - flank)
        end = min(len(seq), positions[read_name]["end"] + flank)
        region = seq[start:end]
        if len(region) > max_read_len:
            continue
        kmers = [region[j:j+kmer_size] for j in range(len(region) - kmer_size + 1)]
        if any(kmer in all_kmers for kmer in kmers):
            for l, kmer in enumerate(kmers):
                if kmer in hits[read_name]:
                    matched_x.append(l)
                    matched_y.append(count)
                    matched_points.setdefault(kmer, []).append((l, count))
                else:
                    unmatched_x.append(l)
                    unmatched_y.append(count)
            count += 1
    return matched_points, matched_x, matched_y, unmatched_x, unmatched_y, hit_kmers

fasta_dict, all_kmers, fasta_dict = load_fasta_and_get_kmers(fasta_file, kmer_size)

matched_points_tp, matched_x_tp, matched_y_tp, unmatched_x_tp, unmatched_y_tp, hit_kmers = load_fastq_and_get_kmers(minimatches_file_tp, tp_reads_file)

matched_points, matched_x, matched_y, unmatched_x, unmatched_y, hit_kmers = load_fastq_and_get_kmers(minimatches_file, reads_file, hit_kmers)

# Reference hits
unmatched_x_ref, unmatched_y_ref = [], []
matched_x_ref, matched_y_ref = [], []
matched_points_ref = {}  # kmer -> list of (x, y)
ref_order = list(fasta_dict.keys())

for i, ref_name in enumerate(ref_order):
    region = fasta_dict[ref_name]
    kmers = [region[j:j+kmer_size] for j in range(len(region) - kmer_size + 1)]
    for l, kmer in enumerate(kmers):
        if kmer in hit_kmers:
            matched_x_ref.append(l)
            matched_y_ref.append(i)
            matched_points_ref.setdefault(kmer, []).append((l, i))
        elif rc(kmer) in hit_kmers:
            matched_x_ref.append(l)
            matched_y_ref.append(i)
            matched_points_ref.setdefault(rc(kmer), []).append((l, i))
        else:
            unmatched_x_ref.append(l)
            unmatched_y_ref.append(i)

# Plotting
fig = plt.figure(figsize=(24, 8))
gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
ax0 = plt.subplot(gs[0])  # reference
ax1 = plt.subplot(gs[1])  # reads
ax2 = plt.subplot(gs[2])  # reads

# Reference panel
ax0.scatter(unmatched_x_ref, unmatched_y_ref, c="lightgrey", label="Unmatched", s=12, alpha=0.3)
colormap = plt.cm.get_cmap("tab20", len(matched_points_ref))
colormapping = {}
for idx, (kmer, coords) in enumerate(matched_points_ref.items()):
    x, y = zip(*coords)
    ax0.scatter(x, y, s=12, alpha=1, color=colormap(idx))
    colormapping[kmer] = colormap(idx)
ax0.set_title("$k$-mer hits to references")
ax0.set_ylabel("reference")
ax0.set_yticks([])
ax0.set_yticklabels([])
ax0.set_xlim([0, 300])
ax0.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

# Reads panel
ax1.scatter(unmatched_x_tp, unmatched_y_tp, c="lightgrey", s=12, label="Unmatched", alpha=0.3)
colormap = plt.cm.get_cmap("tab20", len(matched_points_tp))
for idx, (kmer, coords) in enumerate(matched_points_tp.items()):
    x, y = zip(*coords)
    ax1.scatter(x, y, s=12, alpha=1, label=kmer, color=colormapping[kmer])

ax1.set_xlabel("$k$-mer index")
ax1.set_ylabel("read")
ax1.set_title("$k$-mer hits for true positive calls in the reads")
ax1.set_yticks([])
ax1.set_yticklabels([])
ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
ax1.set_xlim([0, 300])

# Reads panel
ax2.scatter(unmatched_x, unmatched_y, c="lightgrey", s=12, label="Unmatched", alpha=0.3)
colormap = plt.cm.get_cmap("tab20", len(matched_points))
for idx, (kmer, coords) in enumerate(matched_points.items()):
    x, y = zip(*coords)
    ax2.scatter(x, y, s=12, alpha=1, label=kmer, color=colormapping[kmer])

ax2.set_xlabel("$k$-mer index")
ax2.set_ylabel("read")
ax2.set_title("$k$-mer hits for false positive calls in the reads")
ax2.set_yticks([])
ax2.set_yticklabels([])
ax2.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
ax2.set_xlim([0, 300])
plt.tight_layout()

# Save outputs
plt.savefig(f"{outfile_base}.png", dpi=600)
plt.savefig(f"{outfile_base}.pdf")
print(f"Plots saved to {outfile_base}.png/.pdf")

import json
import os
from typing import List
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib import rcParams
import matplotlib.patches as mpatches
import glob
from tqdm import tqdm
import tempfile

sns.set_context("talk")
sns.set_style("white")

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
rcParams['font.size'] = 16
import random
random.seed(42)
np.random.seed(42)

def parse_AMR_fasta(AMR_fasta):
    with open(AMR_fasta) as i:
        rows = i.read().split(">")[1:]
    amr_genes = set()
    for r in rows:
        amr_genes.add(r.split(";")[0])
    return amr_genes

def rv_read(read):
    rv = list(reversed(read))
    return ["+" + g[1:] if g[0] == "-" else "-" + g[1:] for g in rv]

def parse_json(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)

def wrap_around_find_ref_genes(ref_json, amr_genes):
    AMR_loci = {}
    references = {}
    contexts = {}
    for contig in ref_json:
        contig_length = len(ref_json[contig])
        for i, g in enumerate(ref_json[contig]):
            if g[1:] in amr_genes:
                if g[1:] not in AMR_loci:
                    AMR_loci[g[1:]] = {}
                    references[g[1:]] = {}
                    contexts[g[1:]] = {}
                AMR_loci[g[1:]][f"{contig};{i}"] = []
                # Circular prefix
                prefix_start = (i - 30) % contig_length
                if i < 30:
                    prefix = ref_json[contig][prefix_start:] + ref_json[contig][:i]
                else:
                    prefix = ref_json[contig][i-30:i]
                # Circular suffix
                suffix_end = (i + 31) % contig_length
                if i + 31 > contig_length:
                    suffix = ref_json[contig][i+1:] + ref_json[contig][:suffix_end]
                else:
                    suffix = ref_json[contig][i+1:i+31]
                references[g[1:]][f"{contig};{i}"] = prefix + [ref_json[contig][i]] + suffix
                contexts[g[1:]][f"{contig};{i}"] = {"prefix": len(prefix), "suffix": len(suffix)}

    return AMR_loci, references, contexts

def find_ref_genes(ref_json, amr_genes):
    AMR_loci = {}
    references = {}
    contexts = {}
    for contig in ref_json:
        for i, g in enumerate(ref_json[contig]):
            if g[1:] in amr_genes:
                if g[1:] not in AMR_loci:
                    AMR_loci[g[1:]] = {}
                    references[g[1:]] = {}
                    contexts[g[1:]] = {}
                AMR_loci[g[1:]][f"{contig};{i}"] = []
                prefix = ref_json[contig][i-30:i]
                suffix = ref_json[contig][i+1:i+31]
                #if len(prefix) < 30:
                #    prefix = suffix[-10:] + prefix
                #if len(suffix) < 30:
                #    suffix = suffix + prefix[:10]
                references[g[1:]][f"{contig};{i}"] = prefix + [ref_json[contig][i]] + suffix
                contexts[g[1:]][f"{contig};{i}"] = {"prefix": len(prefix), "suffix": len(suffix)}
    return AMR_loci, references, contexts

def assign_AMR_reads(ref_amr_genes, truth_reads):
    for r in truth_reads:
        read_id, contig, start, end = r.split("_")
        start = int(start)
        end = int(end)
        if start < end:
            orientation = "+"
        else:
            orientation = "-"
        start, end = sorted([start, end])
        for g in ref_amr_genes:
            for allele in ref_amr_genes[g]:
                ref_contig, ref_index = allele.split(";")
                if contig != ref_contig:
                    continue
                if start <= int(ref_index) <= end:
                    ref_amr_genes[g][allele].append(orientation + ";" + read_id)
    for g in ref_amr_genes:
        for allele in ref_amr_genes[g]:
            ref_amr_genes[g][allele] = sorted(ref_amr_genes[g][allele])
    return ref_amr_genes

def get_valid_chars():
    excluded = {0x3E, 0x3D, 0x3C, 0x2D, 0x20, 0x0D, 0x0A}
    return [chr(i) for i in range(1, 256) if i not in excluded]

# Updated function to decode alignment characters
def parse_alignment(path: str, char_gene_map):
    with open(path, 'r') as f:
        aln = f.read()
    names = []
    seqs = []
    for entry in aln.split(">")[1:]:
        lines = entry.strip().split("\n")
        name = lines[0]
        hex_seq = " ".join(lines[1:]).split()
        decoded_seq = []
        for hex_char in hex_seq:
            if hex_char == "--":
                decoded_seq.append("-")
            else:
                try:
                    ascii_char = chr(int(hex_char, 16))
                    decoded_seq.append(char_gene_map[ascii_char])
                except (ValueError, KeyError):
                    decoded_seq.append("?")  # fallback for unknown characters
        names.append(name)
        seqs.append(decoded_seq)
    return names, seqs

def make_alignment_matrix(correction_jsons, ref_amr_genes, labels, references, aln_alg, debug_prefix=None):
    valid_reads = set()
    ref_set = set(references)
    unique_genes = set(references)
    for r in ref_amr_genes:
        orientation, read = r.split(";")
        if all(read in j for j in correction_jsons):
            valid_reads.add(r)
    valid_reads = list(sorted(valid_reads))[:200]
    read_data = []
    for i, data in enumerate(correction_jsons):
        label = labels[i].replace(" ", "_")
        for r in valid_reads:
            orientation, read = r.split(";")
            read_label = read + ";" + label
            read_genes = data[read][:]
            rv_read_genes = rv_read(read_genes)
            if len(ref_set.intersection(read_genes)) > len(ref_set.intersection(rv_read_genes)):
                read_data.append((read_label, read_genes))
                unique_genes.update(set(read_genes))
            elif len(ref_set.intersection(read_genes)) < len(ref_set.intersection(rv_read_genes)):
                read_data.append((read_label, rv_read_genes))
                unique_genes.update(set(rv_read_genes))
            else:
                if orientation == "+":
                    read_data.append((read_label, read_genes))
                    unique_genes.update(set(read_genes))
                else:
                    read_data.append((read_label, rv_read_genes))
                    unique_genes.update(set(rv_read_genes))
    valid_chars = get_valid_chars()
    if len(unique_genes) > len(valid_chars):
        raise ValueError("Too many genes for available character space (max 248).")

    gene_char_map = {gene: valid_chars[i] for i, gene in enumerate(unique_genes)}
    char_gene_map = {valid_chars[i]: gene for i, gene in enumerate(unique_genes)}

    def encode_sequence(gene):
        return f"{ord(gene_char_map[gene]):02X}"

    fasta_content = [f">ref\n{' '.join([encode_sequence(g) for g in references])}"]
    for r, genes in read_data:
        fasta_content.append(f">{r}\n{' '.join([encode_sequence(g) for g in genes])}")

    if debug_prefix:
        unaligned_path = debug_prefix + "_unaligned.hex"
        aligned_path =  debug_prefix + "_aligned.hex"
        ascii_unaligned =  debug_prefix + "_unaligned.ASCII"
        ascii_aligned =  debug_prefix + "_aligned.ASCII"
    else:
        tmpdir = tempfile.TemporaryDirectory()
        unaligned_path = os.path.join(tmpdir.name, "unaligned.hex")
        aligned_path = os.path.join(tmpdir.name, "aligned.hex")
        ascii_unaligned = unaligned_path.replace(".hex", ".ASCII")
        ascii_aligned = aligned_path.replace(".hex", ".ASCII")

    with open(unaligned_path, "w") as out:
        out.write("\n".join(fasta_content))

    subprocess.run(f"software/hex2maffttext {unaligned_path} > {ascii_unaligned}", shell=True, check=True)
    subprocess.run(f"mafft --maxiterate 2000 --{aln_alg} --text --op 3 --thread 1 {ascii_unaligned} > {ascii_aligned}", shell=True, check=True)
    subprocess.run(f"software/maffttext2hex {ascii_aligned} > {aligned_path}", shell=True, check=True)

    names, seqs = parse_alignment(aligned_path, char_gene_map)

    if not debug_prefix:
        tmpdir.cleanup()

    ref_seq = seqs[0]
    seqs = seqs[1:]
    names = names[1:]

    no_correction = []
    partial_filtering = []
    node_filtering = []
    iterative_correction = []
    read_counts = {}

    for i, name in enumerate(names):
        base_name = name.split(";")[0]
        if labels[0] in name:
            no_correction.append((base_name, seqs[i]))
        if labels[1] in name:
            partial_filtering.append((base_name, seqs[i]))
        if labels[2] in name:
            node_filtering.append((base_name, seqs[i]))
        if labels[3] in name:
            iterative_correction.append((base_name, seqs[i]))
        if base_name not in read_counts:
            read_counts[base_name] = 0
        read_counts[base_name] += 1

    nc, start_index, end_index = make_matrix(no_correction, ref_seq, read_counts)
    pf, start_index, end_index = make_matrix(partial_filtering, ref_seq, read_counts)
    nf, start_index, end_index = make_matrix(node_filtering, ref_seq, read_counts)
    ic, start_index, end_index = make_matrix(iterative_correction, ref_seq, read_counts)
    return nc, pf, nf, ic, ref_seq, start_index, end_index


def make_matrix(seqs, ref_seq, read_counts):
    start_index = None
    for i, c in enumerate(ref_seq):
        if c != "-" and start_index == None:
            start_index = i
        if c != "-":
            end_index = i
    ref_trimmed = ref_seq[start_index: end_index+1]
    state_mat = np.zeros((len(seqs), len(ref_trimmed)), dtype=int)
    MATCH, INS, DEL, MIS, GAP = 0, 1, 2, 3, 4
    trimmed_state_mat = []
    for r, row in enumerate(seqs):
        new_seq = row[1][start_index: end_index+1]
        assert len(new_seq) == len(ref_trimmed)
        state_row = []
        for c, base in enumerate(new_seq):
            ref_base = ref_trimmed[c]
            if base == "-" and ref_base == "-":
                state = GAP
            elif base == ref_base:
                state = MATCH
            elif ref_base == "-" and base != "-":
                state = INS
            elif ref_base != "-" and base == "-":
                state = DEL
            else:
                state = MIS
            state_row.append(state)
        trimmed_state_mat.append(state_row)

    return np.array(trimmed_state_mat, dtype=int), start_index, end_index


AMR_fasta = "AMR_alleles_unified.fa"
amr_genes = parse_AMR_fasta(AMR_fasta)
no_correction_files = glob.glob("assessment_results/pandora_outputs_no_correction/*/*.json")
output_dir = "correction_plots"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
titles = [
    "no correction",
    "partial gene filtering",
    "node filtering",
    "iterative correction"
]
colors = ["#228833", "#CCBB44", "#BBBBBB", "#66CCEE", "white"]  # MATCH, INS, DEL, MIS
labels = ["Match", "Insertion", "Deletion", "Mismatch", "Gap"]
cmap = ListedColormap(colors)
highlight_color = "#AA3377"
highlight_colormap = ListedColormap([highlight_color])
samples_of_interest = {"CFT073", "Escherichia_coli_MINF_1A"}
heatmap_data = {}

for nc_file in tqdm(no_correction_files):
    sample = os.path.basename(nc_file).replace(".json", "")
    if sample not in samples_of_interest:
        continue

    outdir = os.path.join(output_dir, "combined_samples")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    REF_JSON = f"assessment_results/truth_gff_jsons/{sample}.json"
    truth_reads = parse_json(f"assessment_results/processing/{sample}.json")
    ref_json = parse_json(REF_JSON)
    ref_amr_genes, references, contexts = find_ref_genes(ref_json, amr_genes)
    ref_amr_genes = assign_AMR_reads(ref_amr_genes, truth_reads)

    pre_filtering = os.path.join("assessment_results/amira_outputs", sample, "gene_calls_with_gene_filtering.json")
    node_filtering = os.path.join("assessment_results/amira_outputs", sample, "mid_correction_gene_calls.json")
    iterative_correction = os.path.join("assessment_results/amira_outputs", sample, "corrected_gene_calls.json")
    for gene in ref_amr_genes:
        if gene != "blaEC":
            continue
        for allele in ref_amr_genes[gene]:
            for aln_alg in ["localpair"]:
                contig, index = allele.split(";")
                context = contexts[gene][allele]
                fig_key = f"{gene}_{allele.replace(';', '_')}_{aln_alg}"

                nc, pf, nf, ic, ref_seq, start_index, end_index = make_alignment_matrix([
                    parse_json(nc_file),
                    parse_json(pre_filtering),
                    parse_json(node_filtering),
                    parse_json(iterative_correction)],
                    ref_amr_genes[gene][allele], [t.replace(" ", "_") for t in titles], references[gene][allele], aln_alg)

                column_index = 0
                ref_index = 0
                for b in ref_seq[start_index: end_index+1]:
                    if b != "-":
                        if ref_index == context["prefix"]:
                            break
                        ref_index += 1
                    column_index += 1
                blaec_index = column_index

                # Collect all plots for this gene+allele+alg combo
                if fig_key not in heatmap_data:
                    heatmap_data[fig_key] = []
                heatmap_data[fig_key] = {
                    "sample": sample,
                    "contig": contig,
                    "index": index,
                    "context": context,
                    "matrices": [nc, pf, nf, ic],
                    "ref_seq": ref_seq,
                    "blaec_index": blaec_index
                }

fig, axes = plt.subplots(4, 2, figsize=(16, 24), sharex=True, sharey=True)
axes = axes.flatten()
plotted = False

for col, plot_info in enumerate(heatmap_data):  # col = 0 (CFT073), 1 (MINF_1A)
    sample = heatmap_data[plot_info]["sample"]
    matrices = heatmap_data[plot_info]["matrices"]
    blaec_index = heatmap_data[plot_info]["blaec_index"]
    contig = heatmap_data[plot_info]["contig"]
    index = heatmap_data[plot_info]["index"]
    context = heatmap_data[plot_info]["context"]
    ref_seq = heatmap_data[plot_info]["ref_seq"]

    for row in range(4):  # 4 correction strategies
        ax_idx = row * 2 + col
        ax = axes[ax_idx]
        try:
            sns.heatmap(matrices[row],
                        cmap=cmap,
                        cbar=False,
                        xticklabels=False,
                        yticklabels=False,
                        vmin=0, vmax=4,
                        ax=ax)
            highlight_mask = np.ones_like(matrices[row], dtype=bool)
            highlight_mask[:, blaec_index] = False
            sns.heatmap(np.full_like(matrices[row], 0),
                        cmap=highlight_colormap,
                        mask=highlight_mask,
                        cbar=False,
                        xticklabels=False,
                        yticklabels=False,
                        vmin=0, vmax=1,
                        ax=ax)
            ax.set_title(f"{titles[row]}\n{sample}", fontsize=12)
            plotted = True
        except Exception as e:
            print(f"Skipping subplot {ax_idx} due to error: {e}")
            continue

if plotted:
    fig.supxlabel("Alignment Position", fontsize=16)
    fig.supylabel("Read ID", fontsize=16)

    patches = [mpatches.Patch(facecolor=colors[i], label=labels[i], edgecolor='black', linewidth=0.5) for i in range(5)]
    patches.append(mpatches.Patch(facecolor=highlight_color,
                                    label="$\\mathit{blaEC}$",
                                    edgecolor='black',
                                    linewidth=0.5))
    fig.legend(handles=patches,
                loc='center right',
                frameon=False,
                bbox_to_anchor=(1.15, 0.5))

    plt.tight_layout()
    fig_path_base = os.path.join(output_dir, "combined_samples", f"{fig_key}_combined")
    plt.savefig(f"{fig_path_base}.png", dpi=900, bbox_inches='tight')
    plt.savefig(f"{fig_path_base}.pdf", bbox_inches='tight')
    plt.close(fig)

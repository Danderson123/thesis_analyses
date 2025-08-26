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

sns.set_context("talk")
sns.set_style("white")

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
rcParams['font.size'] = 16


def rv_unitig(gene_list: List[str]) -> List[str]:
    rv = list(reversed(gene_list))
    return ["+" + g[1:] if g[0] == "-" else "-" + g[1:] for g in rv]

def parse_unitigs(path: str, ref_set: set, orientation, central_gene) -> List[List[str]]:
    unitigs = []
    with open(path) as fh:
        for line in fh:
            if central_gene not in line:
                continue
            genes = line.split("\t")[0].split(",")
            rv_genes = rv_unitig(genes)
            if orientation in genes:
                unitigs.append(genes)
            else:
                unitigs.append(rv_genes)
    return unitigs

def rv_read(read):
    rv = list(reversed(read))
    return ["+" + g[1:] if g[0] == "-" else "-" + g[1:] for g in rv]

def parse_json(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)

def make_alignment_matrix(read_files, ref, labels, central_gene):
    unique_genes = set(ref)
    # Load reads
    read_list = []
    read_counts = {}
    for i, read_file in enumerate(read_files):
        read_json = parse_json(read_file)
        for read in read_json:
            if f"+{central_gene}" in read_json[read]:
                if read not in read_counts:
                    read_counts[read] = 0
                read_counts[read] += 1
                read_list.append((read + ";" + labels[i], read_json[read]))
                unique_genes.update(set(read_json[read]))
            elif f"-{central_gene}" in read_json[read]:
                if read not in read_counts:
                    read_counts[read] = 0
                read_counts[read] += 1
                rr = rv_read(read_json[read])
                read_list.append((read + ";" + labels[i], rr))
                unique_genes.update(set(rr))
    import hashlib
    valid_reads = set([k for k, v in read_counts.items() if v == 4][:250])
    #assert len(valid_reads) == 250
    # Step 1: Get the list of valid characters
    def get_valid_chars():
        excluded = {0x3E, 0x3D, 0x3C, 0x2D, 0x20, 0x0D, 0x0A}
        return [chr(i) for i in range(1, 256) if i not in excluded]

    valid_chars = get_valid_chars()
    if len(unique_genes) > len(valid_chars):
        raise ValueError("Too many genes for available character space (max 248).")

    gene_char_map = {gene: valid_chars[i] for i, gene in enumerate(unique_genes)}
    char_gene_map = {valid_chars[i]: gene for i, gene in enumerate(unique_genes)}

    # Step 4: Encode sequences into hex strings
    def encode_sequence(g):
        return f"{ord(gene_char_map[g]):02X}"
    read_list = sorted(read_list)
    # Step 5: Write to input.hex
    fasta_content = [f">ref\n{' '.join([encode_sequence(g) for g in ref])}"]
    for r in read_list:
        if r[0].split(";")[0] in valid_reads:
            fasta_content.append(f">{r[0]}\n{' '.join([encode_sequence(g) for g in r[1]])}")

    unaligned_path = read_file.replace(".json", ".unaligned.hex")
    aligned_path = read_file.replace(".json", ".aligned.hex")
    with open(unaligned_path, "w") as out:
        out.write("\n".join(fasta_content))

    # Convert to ASCII
    command = f"software/hex2maffttext {unaligned_path} > {unaligned_path.replace('.hex', '.ASCII')}"
    subprocess.run(command, shell=True, check=True)

    # Run MAFFT alignment
    command = f'mafft --localpair --text {unaligned_path.replace(".hex", ".ASCII")} > {aligned_path.replace(".hex", ".ASCII")}'
    subprocess.run(command, shell=True, check=True)

    # Convert back to hex
    command = f"software/maffttext2hex {aligned_path.replace('.hex', '.ASCII')} > {aligned_path}"
    subprocess.run(command, shell=True, check=True)

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

    names, seqs = parse_alignment(aligned_path, char_gene_map)
    ref_seq = seqs[0]                              # first row is reference
    seqs = seqs[1:]
    names = names[1:]
    no_correction = []
    partial_filtering = []
    node_filtering = []
    iterative_correction = []
    read_counts = {}
    for i, name in enumerate(names):
        if labels[0] in name:
            name = name.split(";")[0]
            no_correction.append((name, seqs[i]))
        if labels[1] in name:
            name = name.split(";")[0]
            partial_filtering.append((name, seqs[i]))
        if labels[2] in name:
            name = name.split(";")[0]
            node_filtering.append((name, seqs[i]))
        if labels[3] in name:
            name = name.split(";")[0]
            iterative_correction.append((name, seqs[i]))
        if name not in read_counts:
            read_counts[name] = 0
        read_counts[name] += 1
    nc, ref_trimmed= make_matrix(no_correction, ref_seq, read_counts)
    pf, ref_trimmed = make_matrix(partial_filtering, ref_seq, read_counts)
    nf, ref_trimmed = make_matrix(node_filtering, ref_seq, read_counts)
    ic, ref_trimmed = make_matrix(iterative_correction, ref_seq, read_counts)
    return nc, pf, nf, ic, ref_trimmed

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

    return np.array(trimmed_state_mat, dtype=int), ref_trimmed

# Group files by strain and correction stage
cft_files = [
    "assessment_results/pandora_outputs_no_correction/Escherichia_coli_MINF_1A/Escherichia_coli_MINF_1A.json",
    "assessment_results/amira_outputs/Escherichia_coli_MINF_1A/gene_calls_with_gene_filtering.json",
    "assessment_results/amira_outputs/Escherichia_coli_MINF_1A/mid_correction_gene_calls.json",
    "assessment_results/amira_outputs/Escherichia_coli_MINF_1A/corrected_gene_calls.json"
]

minf1a_files = [
    "assessment_results/pandora_outputs_no_correction/CFT073/CFT073.json",
    "assessment_results/amira_outputs/CFT073/gene_calls_with_gene_filtering.json",
    "assessment_results/amira_outputs/CFT073/mid_correction_gene_calls.json",
    "assessment_results/amira_outputs/CFT073/corrected_gene_calls.json"
]

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

fig, axes = plt.subplots(4, 2, figsize=(16, 24), sharex=True, sharey=True)

for strain_label in ["CFT073", "Escherichia_coli_MINF_1A"]:
    REF_JSON = f"assessment_results/truth_gff_jsons/{strain_label}.json"
    # Load reference
    if strain_label == "Escherichia_coli_MINF_1A":
        central_gene = "blaEC"
    else:
        #central_gene = "ant3Ihaac6IIdNG_0472181"
        central_gene = "blaEC"
    ref_json = parse_json(REF_JSON)
    for contig in ref_json:
        if f"+{central_gene}" in ref_json[contig]:
            ref = ref_json[contig]
            orientation = f"+{central_gene}"
        if f"-{central_gene}" in ref_json[contig]:
            ref = rv_read(ref_json[contig])
            orientation = f"-{central_gene}"

    blaec = [i for i, val in enumerate(ref) if f"{central_gene}" in val][0]
    if strain_label == "CFT073":
        ref = ref[blaec - 30: blaec + 31]
        strain_files = cft_files
        i = 0
    else:
        ref = ref[blaec - 30: blaec + 31]
        strain_files = minf1a_files
        i = 1
    nc, pf, nf, ic, ref_seq = make_alignment_matrix(strain_files, ref, [t.replace(" ", "_") for t in titles], central_gene)
    blaec_index = next(k for k, g in enumerate(ref_seq) if f"{central_gene}" in g)
    # plot the heatmaps
    ax = axes[0, i]
    sns.heatmap(nc,
                cmap=cmap,
                cbar=False,
                xticklabels=False,
                yticklabels=False,
                vmin=0, vmax=4,
                ax=ax)
    highlight_mask = np.ones_like(nc, dtype=bool)
    highlight_mask[:, blaec_index] = False
    sns.heatmap(np.full_like(nc, 0),
                cmap=highlight_colormap,
                mask=highlight_mask,
                cbar=False,
                xticklabels=False,
                yticklabels=False,
                vmin=0, vmax=1,
                ax=ax)
    ax.set_title(f"{strain_label} - {titles[0]}")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax = axes[1, i]
    sns.heatmap(pf,
                cmap=cmap,
                cbar=False,
                xticklabels=False,
                yticklabels=False,
                vmin=0, vmax=4,
                ax=ax)
    highlight_mask = np.ones_like(pf, dtype=bool)
    highlight_mask[:, blaec_index] = False
    sns.heatmap(np.full_like(pf, 0),
                cmap=highlight_colormap,
                mask=highlight_mask,
                cbar=False,
                xticklabels=False,
                yticklabels=False,
                vmin=0, vmax=1,
                ax=ax)
    ax.set_title(f"{strain_label} - {titles[1]}")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax = axes[2, i]
    sns.heatmap(nf,
                cmap=cmap,
                cbar=False,
                xticklabels=False,
                yticklabels=False,
                vmin=0, vmax=4,
                ax=ax)
    highlight_mask = np.ones_like(nf, dtype=bool)
    highlight_mask[:, blaec_index] = False
    sns.heatmap(np.full_like(nf, 0),
                cmap=highlight_colormap,
                mask=highlight_mask,
                cbar=False,
                xticklabels=False,
                yticklabels=False,
                vmin=0, vmax=1,
                ax=ax)
    ax.set_title(f"{strain_label} - {titles[2]}")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax = axes[3, i]
    sns.heatmap(ic,
                cmap=cmap,
                cbar=False,
                xticklabels=False,
                yticklabels=False,
                vmin=0, vmax=4,
                ax=ax)
    highlight_mask = np.ones_like(ic, dtype=bool)
    highlight_mask[:, blaec_index] = False
    sns.heatmap(np.full_like(ic, 0),
                cmap=highlight_colormap,
                mask=highlight_mask,
                cbar=False,
                xticklabels=False,
                yticklabels=False,
                vmin=0, vmax=1,
                ax=ax)
    ax.set_title(f"{strain_label} - {titles[3]}")
    ax.set_xlabel("")
    ax.set_ylabel("")

# Shared axis labels
fig.supxlabel("Alignment Position", fontsize=16)
fig.supylabel("Read ID", fontsize=16)

# Legend
patches = [mpatches.Patch(facecolor=colors[i], label=labels[i], edgecolor='black', linewidth=0.5) for i in range(5)]
patches.append(mpatches.Patch(facecolor=highlight_color, label="$\\mathit{blaEC}$", edgecolor='black', linewidth=0.5))

fig.legend(handles=patches,
           loc='center right',
           frameon=False,
           bbox_to_anchor=(1.15, 0.5))

plt.tight_layout()
plt.savefig(os.path.join("assessment_results", "combined_heatmap.bad.png"), dpi=900, bbox_inches='tight')
plt.savefig(os.path.join("assessment_results", "combined_heatmap.bad.pdf"), bbox_inches='tight')


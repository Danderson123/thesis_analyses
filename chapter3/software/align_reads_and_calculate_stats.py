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

def parse_unitigs(path: str, ref_set: set, orientation) -> List[List[str]]:
    unitigs = []
    with open(path) as fh:
        for line in fh:
            if "blaEC" not in line:
                continue
            genes = line.split("\t")[0].split(",")
            rv_genes = rv_unitig(genes)
            if orientation in genes:
                unitigs.append(genes)
            else:
                unitigs.append(rv_genes)
    return unitigs

def parse_json(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)

def make_alignment_matrix(read_files, ref, labels):
    unique_genes = set(ref)
    # Load reads
    read_list = []
    read_counts = {}
    for i, read_file in enumerate(read_files):
        read_json = parse_json(read_file)
        for read in read_json:
            if "+blaEC" in read_json[read]:
                if read not in read_counts:
                    read_counts[read] = 0
                read_counts[read] += 1
                read_list.append((read + ";" + labels[i], read_json[read]))
                unique_genes.update(set(read_json[read]))
            elif "-blaEC" in read_json[read]:
                if read not in read_counts:
                    read_counts[read] = 0
                read_counts[read] += 1
                rr = rv_read(read_json[read])
                read_list.append((read + ";" + labels[i], rr))
                unique_genes.update(set(rr))
    import hashlib
    valid_reads = set([k for k, v in read_counts.items() if v == 4][:50])
    print(len(valid_reads))
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
    command = f'mafft --text {unaligned_path.replace(".hex", ".ASCII")} > {aligned_path.replace(".hex", ".ASCII")}'
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
    nc = make_matrix(no_correction, ref_seq, read_counts)
    pf = make_matrix(partial_filtering, ref_seq, read_counts)
    nf = make_matrix(node_filtering, ref_seq, read_counts)
    ic = make_matrix(iterative_correction, ref_seq, read_counts)
    return nc, pf, nf, ic, ref_seq

def make_matrix(seqs, ref_seq, read_counts):
    state_mat = np.zeros((len(seqs), len(ref_seq)), dtype=int)
    MATCH, INS, DEL, MIS = 0, 1, 2, 3
    trimmed_state_mat = []
    for r, row in enumerate(seqs):
        state_row = []
        for c, base in enumerate(row[1]):
            ref_base = ref_seq[c]
            if base == ref_base:
                state = MATCH
            elif ref_base == "-" and base != "-":
                state = INS
            elif ref_base != "-" and base == "-":
                state = DEL
            else:
                state = MIS
            state_row.append(state)
        trimmed_state_mat.append(state_row)

    return np.array(trimmed_state_mat, dtype=int)

# Group files by strain and correction stage
cft_files = [
    "assessment_results/processing/CFT073.json",
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

def rv_read(read):
    rv = list(reversed(read))
    return ["+" + g[1:] if g[0] == "-" else "-" + g[1:] for g in rv]

def get_valid_chars():
    excluded = {0x3E, 0x3D, 0x3C, 0x2D, 0x20, 0x0D, 0x0A}
    return [chr(i) for i in range(1, 256) if i not in excluded]

def encode_sequence(g):
    return f"{ord(gene_char_map[g]):02X}"

def align_reads(ref, nc, pf, nf, ic, labels, read):
    ref_genes = set(ref)
    unique_genes = set(ref)
    query_list = []
    for q in [nc, pf, nf, ic]:
        rv_q = rv_read(q)
        if len(ref_genes.intersection(q)) >  len(ref_genes.intersection(rv_q)):
            query_list.append(q)
            unique_genes.update(set(q))
        else:
            query_list.append(rv_q)
            unique_genes.update(set(rv_q))
    valid_chars = get_valid_chars()
    if len(unique_genes) > len(valid_chars):
        raise ValueError("Too many genes for available character space (max 248).")
    gene_char_map = {gene: valid_chars[i] for i, gene in enumerate(unique_genes)}
    char_gene_map = {valid_chars[i]: gene for i, gene in enumerate(unique_genes)}
    # Step 5: Write to input.hex
    fasta_content = [f">ref\n{' '.join([encode_sequence(g) for g in ref])}"]
    for i in range(len(query_list)):
        fasta_content.append(f">{read};{labels[i]}\n{' '.join([encode_sequence(g) for g in query_list[i]])}")
    print(fasta_content)

titles = [
    "no correction",
    "partial gene filtering",
    "node filtering",
    "iterative correction"
]

for strain_label in ["CFT073"]:
    ref = {r.split("_")[0]: v for r, v in parse_json(cft_files[0]).items()}
    nc = parse_json(cft_files[1])
    pf = parse_json(cft_files[2])
    nf = parse_json(cft_files[3])
    ic = parse_json(cft_files[4])
    # get the common reads
    common_reads = set(ref.keys()) & set(nc.keys()) & set(pf.keys()) & set(nf.keys()) & set(ic.keys())
    # iterate through the reads
    for read in common_reads:
        alignment = align_reads(ref[read], nc[read], pf[read], nf[read], ic[read], [t.replace(" ", "_") for t in titles], read)
    ss
    nc, pf, nf, ic, ref_seq = make_alignment_matrix(strain_files, ref, [t.replace(" ", "_") for t in titles])

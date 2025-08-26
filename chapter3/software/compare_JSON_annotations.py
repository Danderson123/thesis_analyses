from itertools import product
from collections import deque
import argparse
import json
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import numpy as np

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Get the Pandora list of genes per real read for a specified GFF.')
    parser.add_argument('--truth-annotation', dest='truth_annotation',
                        help='path to truth annotation JSON', required=True)
    parser.add_argument('--pandora-annotation', dest='pandora_annotation',
                        help='path to Pandora annotation JSON', required=True)
    parser.add_argument('--gene-lengths', dest='gene_lengths',
                        help='path to gene length JSON', required=True)
    parser.add_argument('--output-stats', dest='output_stats',
                        help='path to output file of summary statistics', required=True)
    args = parser.parse_args()
    return args

def needleman_wunsch(x, y):
    N, M = len(x), len(y)

    # Modify the scoring function to handle combined elements or sequences
    def s(a, b):
        if "/" in a or "/" in b:  # Check if elements are combined sequences
            a_parts = a.split("/") if "/" in a else [a, '*']
            b_parts = b.split("/") if "/" in b else [b, '*']
            score = 0
            for a_part, b_part in zip(a_parts, b_parts):
                score += int(a_part == b_part or '*' in [a_part, b_part])  # Consider '*' as a wildcard matching any character
            return score
        else:
            return int(a == b)  # Regular comparison for non-combined elements

    DIAG = -1, -1
    LEFT = -1, 0
    UP = 0, -1

    F = {}
    Ptr = {}
    F[-1, -1] = 0
    for i in range(N):
        F[i, -1] = -i
    for j in range(M):
        F[-1, j] = -j

    for i, j in product(range(N), range(M)):
        option_F = (
            F[i - 1, j - 1] + s(x[i], y[j]),
            F[i - 1, j] - 1,
            F[i, j - 1] - 1,
        )
        F[i, j], Ptr[i, j] = max(zip(option_F, (DIAG, LEFT, UP)))

    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = (x[i], y[j])
        elif direction == LEFT:
            element = (x[i], "*")
        elif direction == UP:
            element = ("*", y[j])
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj

    while i >= 0:
        alignment.appendleft((x[i], "*"))
        i -= 1
    while j >= 0:
        alignment.appendleft(("*", y[j]))
        j -= 1

    return list(alignment)

def align_fast(x, y):
    return needleman_wunsch(x, y)

def load_json(filename):
    with open(filename, 'r') as f:
        content = json.loads(f.read())
    new_content = {}
    read_mapping = {}
    for read in content:
        new_content[read.split("_")[0]] = content[read]
        read_mapping[read.split("_")[0]] = read
    return new_content, read_mapping

def align_genes(file1, file2, output_stats):
    data1, read_mapping = load_json(file1)
    data2, _ = load_json(file2)
    sample = os.path.basename(os.path.dirname(file2))
    data3, _ = load_json(f"/hps/nobackup/iqbal/dander/thesis_figures/correction_assessment/assessment_results/amira_outputs/{sample}/corrected_gene_calls.json")
    common_reads = set(data1.keys()) & set(data2.keys()) & set(data3.keys())
    all_alignments = {}
    alignment_lines = []
    for read in tqdm(common_reads):
        data1[read] = [g.replace(".aln.fas", "") for g in data1[read]]
        data2[read] = [g.replace(".aln.fas", "") for g in data2[read]]
        alignment = align_fast(data1[read], data2[read])
        all_alignments[read] = alignment
        # get the alignment rows as strings
        truth = [col[0] for col in alignment]
        pandora = [col[1] for col in alignment]
        alignment_lines.append(">" + read_mapping[read] + "\n" + ",".join(truth) + "\n" + ",".join(pandora))
    # write out the alignments to a file
    file_prefix = output_stats.replace(".txt", "").replace(".summary_stats", "")
    with open(f"{file_prefix}_alignments.txt", "w") as o:
        o.write("\n\n".join(alignment_lines))
    return all_alignments

def plot_top_missed_genes(gene_dict, filename):
    # Sort the dictionary by values (counts) in descending order and select the top 50
    top_genes = sorted(gene_dict.items(), key=lambda x: x[1], reverse=True)[:50]
    # Unzip into two lists: gene names and counts
    genes, counts = zip(*top_genes)
    # Create a bar plot
    plt.figure(figsize=(10, 8))  # Adjust the size as needed
    plt.bar(genes, counts, color='skyblue')
    plt.xlabel('Gene Name')
    plt.ylabel('Missed Count')
    plt.xticks(rotation=90)  # Rotate gene names for better readability
    plt.title('Top 50 Most Frequently Missed Genes')
    plt.tight_layout()  # Adjusts the plot to ensure everything fits without overlapping
    # Save the plot to a PDF file
    plt.savefig(filename, format='pdf', bbox_inches='tight')

def plot_top_fp_genes(gene_dict, pandora_called_counts, filename):
    # print the mean fp rate
    import statistics
    print(f"FP calls: {statistics.mean(list(gene_dict.values()))}")
    # Sort the dictionary by values (counts) in descending order and select the top 50
    top_genes = sorted(gene_dict.items(), key=lambda x: x[1], reverse=True)[:50]
    # Unzip into two lists: gene names and counts
    genes, counts = zip(*top_genes)
    for i in range(len(genes)):
        print(genes[i], counts[i] * 100 / pandora_called_counts[genes[i]], counts[i], pandora_called_counts[genes[i]])
    # Create a bar plot
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.bar(genes, counts, color='skyblue')
    ax.set_xlabel('Gene name')
    ax.set_ylabel('False positive count')
    ax.set_ylim([0, 100])
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90)
    for label in ax.get_xticklabels():
        label.set_fontstyle('italic')

    # Beautify the axes
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)

    plt.tight_layout()
    plt.savefig(filename, format='pdf', bbox_inches='tight')
    plt.close()

def plot_FNs_against_length(gene_missed_counts, gene_read_calls, gene_lengths, output_pdf):
    # load the length json
    with open(gene_lengths) as i:
        gene_lengths = json.load(i)
    # get the length of each gene
    FN = []
    lengths = []
    for g in gene_missed_counts:
        rate = gene_missed_counts[g] / (gene_read_calls[g]["correct"] + gene_read_calls[g]["incorrect"]) * 100
        if int(rate) > 0:
            FN.append(rate)
            lengths.append(gene_lengths[g])
    # make the plot
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(lengths, FN)
    ax.set_xlabel('Gene length')
    ax.set_ylabel('False negative rate')
    ax.set_xlim([0, 5000])
    ax.set_ylim([0, 1])
    ax.set_xticks(np.arange(0, 5500, 500))
    ax.set_yticks(np.arange(0, 110, 10))
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
    plt.tight_layout()  # Adjusts the plot to ensure everything fits without overlapping
    # Save the plot to a PDF file
    plt.savefig(output_pdf, format='pdf', bbox_inches='tight')
    plt.close()

def plot_combined(gene_dict, gene_missed_counts, gene_read_calls, gene_lengths, output_pdf):
    # Prepare data for the first plot
    top_genes = sorted(gene_dict.items(), key=lambda x: x[1], reverse=True)[:50]
    genes, counts = zip(*top_genes)
    with open(output_pdf.replace(".pdf", ".fp_genes.txt"), "w") as o:
        o.write("\n".join(genes))
    # Prepare data for the second plot
    with open(gene_lengths) as i:
        gene_lengths_data = json.load(i)

    FN = []
    lengths = []
    for g in gene_missed_counts:
        rate = gene_missed_counts[g] / (gene_read_calls[g]["correct"] + gene_read_calls[g]["incorrect"]) * 100
        if int(rate) > 0:
            FN.append(rate)
            lengths.append(gene_lengths_data[g])

    # Create the subplots
    fig, ax = plt.subplots(2, 1, figsize=(22, 22))  # One row, two columns

    # Plot 1: Top falsely called genes (bar plot)
    ax[0].bar(genes, counts, color='skyblue', edgecolor='black')
    ax[0].set_xlabel('Gene name', fontsize=18)
    ax[0].set_ylabel('False positive count', fontsize=18)
    ax[0].tick_params(axis='both', labelsize=18)
    ax[0].spines['left'].set_linewidth(0.5)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['bottom'].set_linewidth(0.5)
    ax[0].grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
    ax[0].tick_params(axis='x', rotation=90)
    for label in ax[0].get_xticklabels():
        label.set_fontstyle('italic')

    # Plot 2: FN rate vs. gene length (scatter plot)
    ax[1].scatter(lengths, FN, color='skyblue', s=56, edgecolor='black')
    ax[1].set_xlabel('Gene length', fontsize=18)
    ax[1].set_ylabel('False negative rate (%)', fontsize=18)
    ax[1].tick_params(axis='both', labelsize=18)
    ax[1].spines['left'].set_linewidth(0.5)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['bottom'].set_linewidth(0.5)
    ax[1].grid(axis="y", linestyle="--", alpha=0.7, zorder=1)
    ax[1].set_xlim([0, 5000])
    ax[1].set_ylim([0, 110])
    ax[1].set_xticks(np.arange(0, 5500, 500))
    ax[1].set_yticks(np.arange(0, 110, 10))

    # Adjust layout
    plt.tight_layout()
    # Save the combined plot to a PDF
    plt.savefig(output_pdf, format='pdf', bbox_inches='tight')
    plt.close()

def plot_top_incorrectly_called_genes(gene_dict, filename, gene_read_calls):
    # Sort the dictionary by values (counts) in descending order and select the top 50
    top_genes = sorted(gene_dict.items(), key=lambda x: x[1], reverse=True)[:50]
    # Unzip into two lists: gene names and counts
    genes, counts = zip(*top_genes)
    # get the name of the most incorrectly called gene
    for i in [0, 1, 2]:
        most_incorrect = gene_read_calls[genes[i]]
        # write the list of reads for the correct and incorrect calls
        with open(os.path.join(filename.replace("_most_incorrectly_called_genes.pdf", "") + f"_{genes[i]}_correct_reads.txt"), "w") as o:
            o.write("\n".join(most_incorrect["correct"]))
        with open(os.path.join(filename.replace("_most_incorrectly_called_genes.pdf", "") + f"_{genes[i]}_incorrect_reads.txt"), "w") as o:
            o.write("\n".join(most_incorrect["incorrect"]))
    # Create a bar plot
    plt.figure(figsize=(10, 8))  # Adjust the size as needed
    plt.bar(genes, counts, color='skyblue')
    plt.xlabel('Gene Name')
    plt.ylabel('Count')
    plt.xticks(rotation=90)  # Rotate gene names for better readability
    plt.title('Top 50 Most Frequently Incorrectly Called Genes')
    plt.tight_layout()  # Adjusts the plot to ensure everything fits without overlapping
    # Save the plot to a PDF file
    plt.savefig(filename, format='pdf', bbox_inches='tight')

def gather_alignment_stats(alignment, output_stats):
    # initialise the stat counts
    match_count = 0
    snp_error = 0
    strand_error = 0
    insertion_error = 0
    deletion_error = 0
    insertion_at_end = 0
    deletion_at_end = 0
    total_count = 0
    # initialise a dictioary to keep track of how often each individual gene is missed
    gene_DEL_counts = {}
    # initialise a dictioary to keep track of how often each individual gene is falsely called
    gene_IN_counts = {}
    # initialise a dictionary to keep track of how often each individual gene is called another gene
    gene_SNP_counts = {}
    # initialise a dictioary to keep track of the reads were each gene is correctly or incorrectly called
    gene_call_counts = {}
    # initialise a dictionary to keep track of the number of times each gene is called
    pandora_called_counts = {}
    # count the number and type of errors we see in the alignment
    false_negative_counts = {}
    group_2237 = set()
    insB = set()
    for read in tqdm(alignment):
        called_non_gaps = [i for i, val in enumerate(alignment[read]) if alignment[read][i][1] != "*"]
        call_start = max(1, called_non_gaps[0])
        call_end = min(called_non_gaps[-1], len(alignment[read]) - 2)
        for index, col in enumerate(alignment[read]):
            match = False
            if col[0][1:] != "" and col[0][1:] not in gene_call_counts:
                gene_call_counts[col[0][1:]] = {"correct": 0, "incorrect": 0}
                gene_DEL_counts[col[0][1:]] = 0
                gene_SNP_counts[col[0][1:]] = 0
            if col[1][1:] != "" and col[1][1:] not in gene_IN_counts:
                gene_IN_counts[col[1][1:]] = 0
            # if no indel error
            if col[0] != "*" and col[1] != "*":
                if col[0] != col[1]:
                    # if a snp error
                    if col[0][1:] != col[1][1:]:
                        snp_error += 1
                        gene_SNP_counts[col[0][1:]] += 1
                    # if a strand error
                    else:
                        strand_error += 1
                # if a match
                else:
                    match_count += 1
                    gene_call_counts[col[0][1:]]["correct"] += 1
                    match = True
            else:
                # if not at read ends
                if index >= call_start and index <= call_end:
                    if col[0] != "*" and col[1] == "*":
                        deletion_error += 1
                        gene_DEL_counts[col[0][1:]] += 1
                    if col[0] == "*" and col[1] != "*":
                        insertion_error += 1
                        gene_IN_counts[col[1][1:]] += 1
                        if "group_2237" in col[1]:
                            group_2237.add(read)
                        if "insB" in col[1]:
                            insB.add(read)
                # if at read end
                else:
                    if col[0] != "*" and col[1] == "*":
                        deletion_at_end += 1
                    if col[0] == "*" and col[1] != "*":
                        insertion_at_end += 1
            if col[0] != "*" and match is not True:
                gene_call_counts[col[0][1:]]["incorrect"] += 1
            # track how often pandora calls each gene
            if col[1] != "*":
                if col[1][1:] not in pandora_called_counts:
                    pandora_called_counts[col[1][1:]] = 0
                pandora_called_counts[col[1][1:]] += 1
    total_count = sum([match_count, snp_error, strand_error, insertion_error, deletion_error, insertion_at_end, deletion_at_end])
    return match_count, snp_error, strand_error, insertion_error, deletion_error, insertion_at_end, deletion_at_end, total_count, gene_call_counts, gene_DEL_counts, gene_SNP_counts, gene_IN_counts, pandora_called_counts

if __name__=="__main__":
    # get command line args
    args = get_options()
    # get the alignment of the reads shared between the two files
    alignment = align_genes(args.truth_annotation, args.pandora_annotation, args.output_stats)
    # gather the alignment stats
    match_count, snp_error, strand_error, insertion_error, deletion_error, insertion_at_end, deletion_at_end, total_count, gene_call_counts, gene_DEL_counts, gene_SNP_counts, gene_IN_counts, pandora_called_counts = gather_alignment_stats(alignment, args.output_stats)
    # plot the counts of the top 50 most missed genes
    file_prefix = args.output_stats.replace(".txt", "").replace(".summary_stats", "")
    # plot gene length against FNs
    plot_FNs_against_length(gene_DEL_counts, gene_call_counts, args.gene_lengths, f"{file_prefix}_FNs_against_length.pdf")
    # plot the most frequently called false positive genes
    plot_top_fp_genes(gene_IN_counts, pandora_called_counts, f"{file_prefix}_FP_genes.pdf")
    # combine the plots above
    plot_combined(gene_IN_counts, gene_DEL_counts, gene_call_counts, args.gene_lengths, f"{file_prefix}_combined_plot.pdf")
    # write out the alignment stats
    lines = [
        f"Match: {match_count / total_count * 100}",
        f"SNP error: {snp_error / total_count * 100}",
        f"Deletion error: {deletion_error /  total_count * 100}",
        f"Insertion error: {insertion_error /  total_count * 100}",
        f"Insertion error at end: {insertion_at_end /  total_count * 100}",
        f"Deletion error at end: {deletion_at_end /  total_count * 100}",
        f"Strand error: {strand_error / total_count * 100}",
        f"Per gene recall: {match_count / (match_count + snp_error + deletion_error + strand_error)}",
        f"Per gene precision: {match_count / (match_count + snp_error + insertion_error + strand_error)}",
        ]
    with open(args.output_stats, "w") as o:
        o.write("\n".join(lines))

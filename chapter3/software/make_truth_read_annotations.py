import argparse
import json
import os
import pysam
import subprocess
from tqdm import tqdm

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Get the true list of genes per real read for a specified GFF.')
    parser.add_argument('--bam-file', dest='bam_file',
                        help='path to output BAM', required=True)
    parser.add_argument('--output-json', dest='output_json',
                        help='path to output JSON of true read annotations', required=True)
    parser.add_argument('--fastq-file', dest='fastq_file',
                        help='path to fastq file of simulated reads', required=True)
    parser.add_argument('--gff-file', dest='gff_file',
                        help='path to simulated GFF file', required=True)
    parser.add_argument('--threads', dest='threads',
                        help='Number of CPUs', type=int, required=True)
    parser.add_argument('--panaroo-gff', dest='panaroo_gff', action='store_true')
    parser.add_argument('--bakta-gff', dest='bakta_gff', action='store_true')
    args = parser.parse_args()
    return args

def parse_gff(gff_file, panaroo_gff, bakta_gff):
    '''Parse the GFF file and return a dictionary with keys being sequence ids and values being a list of genes.'''
    genes = {}
    no_name_id = 0
    with open(gff_file, 'r') as gff:
        for line in tqdm(gff):
            if line.startswith("##FASTA"):
                break
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
                if feature_type == "CDS" or feature_type == "Panaroo CDS":
                    # Get the gene ID from the attributes
                    gene_name = attributes.split(";")[0].replace("Name=", "")
                    if seqid not in genes:
                        genes[seqid] = []
                    genes[seqid].append((int(start), int(end), strand + gene_name))
    return genes

def get_genes_from_read(genes_dict, read_id, mapping_start, mapping_end):
    '''Return list of genes for the given range in the given read.'''
    found_genes = []
    found_genes_with_indices = []
    start_index = 10000000
    end_index = -1
    if not read_id in genes_dict:
        return None, None, None
    for i in range(len(genes_dict[read_id])):
        gene_start = genes_dict[read_id][i][0]
        gene_end = genes_dict[read_id][i][1]
        gene_id = genes_dict[read_id][i][2]
        if (mapping_end >= gene_end and mapping_start <= gene_start):
            found_genes.append(gene_id)
            if i < start_index:
                start_index = i
            if i > end_index:
                end_index = i
    return found_genes, start_index, start_index + len(found_genes) - 1

def reverse_gene_list(genes):
    reversed_genes = list(reversed(genes))
    for i in range(len(reversed_genes)):
        if reversed_genes[i][0] == "-":
            reversed_genes[i] = "+" + reversed_genes[i][1:]
        else:
            reversed_genes[i] = "-" + reversed_genes[i][1:]
    return reversed_genes

def get_genes_on_read(gff_dict, read, ref_start, ref_end):
    # get the genes that are entirely contained within the read
    genes, start_index, end_index = get_genes_from_read(gff_dict,
                                                    read.reference_name,
                                                    ref_start,
                                                    ref_end)
    if genes:
        # reverse the list of genes if the read has mapped in the reverse direction
        if read.is_reverse:
            genes = reverse_gene_list(genes)
            read_start = end_index
            read_end = start_index
        else:
            read_start = start_index
            read_end = end_index
        return genes, read_start, read_end
    else:
        return None, None, None

if __name__=="__main__":
    # parse command line arguements
    args = get_options()
    # map the simulated reads to the assembly
    sam_file = args.bam_file.replace(".bam", ".sam")
    subprocess.run(f"minimap2 -x map-ont --eqx --MD -a -t {args.threads} -o {sam_file} {args.gff_file} {args.fastq_file}", shell=True, check=True)
    subprocess.run(f"samtools sort -t {args.threads} -o {args.bam_file} {sam_file}", shell=True, check=True)
    subprocess.run(f"samtools index -@ {args.threads} {args.bam_file}", shell=True, check=True)
    # parse the gff file
    gff_dict = parse_gff(args.gff_file, args.panaroo_gff, args.bakta_gff)
    # load the BAM
    try:
        bamfile = pysam.AlignmentFile(args.bam_file, "rb")
    except:
        os.remove(args.bam_file)
        os.remove(sam_file)
        os.remove(args.bam_file + ".bai")
        sys.exit(1)
    # Map reads to genes
    read_gene_map = {}
    readcounts = {}
    reads = [read for read in tqdm(bamfile.fetch())]
    kept_reads = set()
    for read in tqdm(reads):
        read_id = read.query_name
        # skip getting the genes if the read is unmapped or this is not a primary mapping
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        # get the proportion of matching bases
        matching_bases = sum(length for op, length in read.cigartuples if op == 7)  # op == 0 corresponds to 'M' (match/mismatch)
        total_bases = read.query_length  # Total length of the read
        match_proportion = matching_bases / total_bases if total_bases > 0 else 0
        if match_proportion < 0.9:
            continue
        # check if the read fully maps within the reference
        genes, read_start, read_end = get_genes_on_read(gff_dict, read, read.reference_start, read.reference_end)
        if genes:
            # store the genes
            read_gene_map[f"{read_id}_{read.reference_name}_{str(read_start)}_{read_end}"] = genes[:]
    bamfile.close()
    # write out the true gene annotations
    with open(args.output_json, "w") as o:
        o.write(json.dumps(read_gene_map))

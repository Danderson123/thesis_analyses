import argparse
import glob
import gzip
import json
import os
import pysam

def get_options():
    """define args from the command line"""
    parser = argparse.ArgumentParser(description='Get the Pandora list of genes per real read for a specified GFF.')
    parser.add_argument('--pandora-dir', dest='pandora_dir',
                        help='path to Pandora output', required=True)
    parser.add_argument('--output-json', dest='output_json',
                        help='path to output JSON of Pandora read annotations', required=True)
    args = parser.parse_args()
    return args

def determine_gene_strand(read):
    strandlessGene = read.reference_name.replace(".aln.fas", "").replace(".fasta", "").replace(".fa", "")#.replace("~~~", ";")
    if not read.is_forward:
        gene_name = "-" + strandlessGene
    else:
        gene_name = "+" + strandlessGene
    return gene_name, strandlessGene

def parse_fastq_lines(fh):
    # Initialize a counter to keep track of the current line number
    line_number = 0
    # Iterate over the lines in the file
    for line in fh:
        # Increment the line number
        line_number += 1
        # If the line number is divisible by 4, it's a sequence identifier line
        if line_number % 4 == 1:
            # Extract the identifier from the line
            identifier = line.split(" ")[0][1:]
        # If the line number is divisible by 4, it's a sequence line
        elif line_number % 4 == 2:
            sequence = line.strip()
        elif line_number % 4 == 0:
            # Yield the identifier, sequence and quality
            yield identifier, sequence, line.strip()

def parse_fastq(fastq_file):
    # Initialize an empty dictionary to store the results
    results = {}
    # Open the fastq file
    if ".gz" in fastq_file:
        with gzip.open(fastq_file, 'rt') as fh:
            # Iterate over the lines in the file
            for identifier, sequence, quality in parse_fastq_lines(fh):
                # Add the identifier and sequence to the results dictionary
                results[identifier] = {"sequence": sequence,
                                    "quality": quality}
    else:
        with open(fastq_file, 'r') as fh:
            # Iterate over the lines in the file
            for identifier, sequence, quality in parse_fastq_lines(fh):
                # Add the identifier and sequence to the results dictionary
                results[identifier] = {"sequence": sequence,
                                    "quality": quality}
    # Return the dictionary of results
    return results

def convert_pandora_output(pandoraSam,
                        pandora_consensus):
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(pandoraSam, "rb")
    annotatedReads = {}
    # iterate through the read regions
    for read in pandora_sam_content.fetch():
        # convert the cigarsting to a Cigar object
        cigar = read.cigartuples
        # check if the read has mapped to any regions
        if read.is_mapped:
            # append the strand of the match to the name of the gene
            gene_name, strandlessGene = determine_gene_strand(read)
            # exclude genes that do not have a pandora consensus
            if strandlessGene in pandora_consensus:
                if not read.query_name in annotatedReads:
                    annotatedReads[read.query_name] = []
            # store the per read gene names
            if not read.query_name in annotatedReads:
                annotatedReads[read.query_name] = []
            annotatedReads[read.query_name].append(gene_name)
    assert not len(annotatedReads) == 0
    return annotatedReads

if __name__=="__main__":
    # get command line arguements
    args = get_options()
    # parse the pandora consensus
    pandora_consensus = parse_fastq(glob.glob(os.path.join(args.pandora_dir, "*.fq.gz"))[0])
    # convert the pandora sam to a dictionary
    genes_on_reads = convert_pandora_output(glob.glob(os.path.join(args.pandora_dir, "*.sam"))[0],
                                        pandora_consensus)
    # write out the genes on reads as a json
    with open(args.output_json, "w") as o:
        o.write(json.dumps(genes_on_reads))
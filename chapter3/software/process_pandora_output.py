import glob
import os
from tqdm import tqdm
import json
import pysam

def determine_gene_strand(read: pysam.libcalignedsegment.AlignedSegment) -> tuple[str, str]:
    strandlessGene = (
        read.reference_name.replace("~~~", ";")
        .replace(".aln.fas", "")
        .replace(".fasta", "")
        .replace(".fa", "")
    )
    if not read.is_forward:
        gene_name = "-" + strandlessGene
    else:
        gene_name = "+" + strandlessGene
    return gene_name, strandlessGene

def convert_pandora_output(
    pandoraSam: str,
) -> tuple[dict[str, list[str]], list[str]]:
    # load the pseudo SAM
    pandora_sam_content = pysam.AlignmentFile(pandoraSam, "rb")
    annotatedReads: dict[str, list[str]] = {}
    gene_position_dict: dict[str, list[tuple[int, int]]] = {}
    geneCounts: dict[str, int] = {}
    for read in pandora_sam_content.fetch():
        # convert the cigarsting to a Cigar object
        cigar = read.cigartuples
        # check if the read has mapped to any regions
        if read.is_mapped:
            # get the start base that the region maps to on the read
            # append the strand of the match to the name of the gene
            gene_name, strandlessGene = determine_gene_strand(read)
            # exclude genes that do not have a pandora consensus
            read_name = read.query_name
            if read_name not in annotatedReads:
                annotatedReads[read_name] = []
            # count how many times we see each gene
            if strandlessGene not in geneCounts:
                geneCounts[strandlessGene] = 0
            geneCounts[strandlessGene] += 1
            # store the per read gene names
            annotatedReads[read_name].append(gene_name)
    assert not len(annotatedReads) == 0
    return annotatedReads

outdir = "assessment_results/pandora_outputs_no_correction"
if not os.path.exists(outdir):
    os.mkdir(outdir)
for f in tqdm(glob.glob("assessment_results/amira_outputs/*/pandora_output/*.nanopore.fastq.filtered.sam")):
    sample = os.path.basename(os.path.basename(f)).replace(".nanopore.fastq.filtered.sam", "")
    annotatedReads = convert_pandora_output(f)
    if not os.path.exists(os.path.join(outdir, sample)):
        os.mkdir(os.path.join(outdir, sample))
    with open(os.path.join(outdir, sample, f"{sample}.json"), "w") as o:
        o.write(json.dumps(annotatedReads))
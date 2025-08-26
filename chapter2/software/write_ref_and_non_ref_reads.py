import os
import glob
import pysam

def get_longest_contig(bamfile):
    """Return the name of the longest contig in the BAM header."""
    references = bamfile.header['SQ']
    longest = max(references, key=lambda r: r['LN'])
    return longest['SN']

def process_bam_file(bam_path):
    """Process a BAM file and return two sets of read names: longest contig and other contigs."""
    bam = pysam.AlignmentFile(bam_path, "rb")

    longest_contig = get_longest_contig(bam)

    reads_on_longest = set()
    reads_on_others = set()

    for read in bam:
        if read.is_unmapped:
            continue
        ref_name = bam.get_reference_name(read.reference_id)
        if ref_name == longest_contig:
            reads_on_longest.add(read.query_name)
        else:
            reads_on_others.add(read.query_name)

    bam.close()
    return reads_on_longest, reads_on_others, longest_contig

def write_read_list(reads, filename):
    """Write read names to a file, one per line."""
    with open(filename, 'w') as f:
        for read in sorted(reads):
            f.write(f"{read}\n")

def main(directory):
    bam_files = glob.glob(os.path.join(directory, "*.bam"))
    chromosome_reads = set()
    other_contig_reads = set()
    for bam_path in bam_files:
        print(f"Processing {bam_path}...")
        reads_longest, reads_others, longest_contig = process_bam_file(bam_path)
        chromosome_reads.update(reads_longest)
        other_contig_reads.update(reads_others)

    base = os.path.dirname(directory)
    longest_file = f"{base}/reads_longest_contig.txt"
    others_file = f"{base}/reads_other_contigs.txt"

    write_read_list(chromosome_reads, longest_file)
    write_read_list(other_contig_reads, others_file)

if __name__ == "__main__":
    main("pandora_assessment_with_AMR/0.8_0/processing")

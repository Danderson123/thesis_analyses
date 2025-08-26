import os
import subprocess
import pandas as pd
import pysam

def run_minimap2(assembly_file, genes_file, output_file):
    cmd = f"minimap2 -a -x asm20 --MD -t 12 --eqx {assembly_file} {genes_file} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)

def parse_sam(sam_file, similarity_threshold=90.0, length_threshold=95.0):
    hits = []
    samfile = pysam.AlignmentFile(sam_file, "r")
    for read in samfile.fetch():
        if read.is_unmapped or read.is_supplementary or "SUP" in read.query_name:
            continue
        query_name = read.query_name
        query_length = read.infer_read_length()
        query_start = read.query_alignment_start
        query_end = read.query_alignment_end
        strand = '-' if read.is_reverse else '+'
        target_name = read.reference_name
        target_start = read.reference_start
        target_end = read.reference_end
        alignment_block_length = read.query_alignment_length
        mapq = read.mapping_quality
        # Calculate matching bases and the length used in the alignment
        matching_bases = 0
        total_aligned_length = 0
        for op, length in read.cigartuples:
            if op == 7:  # Match/Mismatch
                matching_bases += length
            if op != 1 and op != 4 and op != 5:  # Not an insertion to reference
                total_aligned_length += length
        # Calculate the similarity and the alignment length percentage
        similarity = (matching_bases / total_aligned_length) * 100.0 if total_aligned_length > 0 else 0
        alignment_length_percentage = (total_aligned_length / query_length) * 100.0 if query_length > 0 else 0
        # Check against thresholds
        if similarity >= similarity_threshold and alignment_length_percentage >= length_threshold:
            hits.append({
                'query_name': query_name,
                'query_length': query_length,
                'query_start': query_start,
                'query_end': query_end,
                'strand': strand,
                'target_name': target_name,
                'target_start': target_start,
                'target_end': target_end,
                'similarity': similarity,
                'matching_bases': matching_bases,
                'alignment_block_length': alignment_block_length,
                'mapq': mapq
            })
    samfile.close()
    return hits

def filter_longest_hits(hits):
    hits_dict = {}
    for row in hits:
        for k, v in row.items():
            if k not in hits_dict:
                hits_dict[k] = []
            hits_dict[k].append(v)
    hits_df = pd.DataFrame(hits_dict)
    hits_df['length'] = hits_df['query_end'] - hits_df['query_start']

    def filter_function(df):
        df = df.sort_values(by=['similarity', 'length'], ascending=[False, False])
        longest_hits = []
        for _, row in df.iterrows():
            overlap = False
            for hit in longest_hits:
                if not (row['target_end'] < hit['target_start'] + 200 or row['target_start'] > hit['target_end'] - 200):
                    overlap = True
                    break
            if not overlap:
                longest_hits.append(row)
        return pd.DataFrame(longest_hits)

    filtered_hits = hits_df.groupby('target_name').apply(filter_function).reset_index(drop=True)
    # Select only one hit per gene type, the longest one with highest similarity if lengths are the same
    return filtered_hits

def write_json(filtered_hits, output_json, assembly_file, additional_genes):
    import json
    filtered_hits = filtered_hits.sort_values(by=['target_name', 'target_start'])
    gene_counts = {}
    for _, hit in filtered_hits.iterrows():
        gene_name = apply_rules(hit['query_name'].split(';')[1].split(".")[0])
        if gene_name not in gene_counts:
            gene_counts[gene_name] = 0
        gene_counts[gene_name] += 1
    for gene in additional_genes:
        gene_name = apply_rules(gene.split(";")[1].split(".")[0])
        if gene_name not in gene_counts:
            gene_counts[gene_name] = 1
    with open(output_json, "w") as o:
        o.write(json.dumps(gene_counts))

def write_gff(filtered_hits, output_file, assembly_file):
    filtered_hits = filtered_hits.sort_values(by=['target_name', 'target_start'])
    with open(assembly_file) as i:
        assembly_content = i.read()
    with open(output_file, 'w') as f:
        for _, hit in filtered_hits.iterrows():
            gff_line = f"{hit['target_name']}\tminimap2\tCDS\t{hit['target_start'] + 1}\t{hit['target_end']}\t.\t{hit['strand']}\t.\tName={';'.join(hit['query_name'].split(';')[:2])}\n"
            f.write(gff_line)
        f.write(f"##FASTA\n{assembly_content}")

def supplement_with_genes_from_reads(genes_file, nanopore, output_file):
    # map the reads to the genes
    cmd = f"minimap2 -a -x map-ont -t 12 --eqx --MD {genes_file} {nanopore} > {output_file}"
    if not os.path.exists(output_file):
        subprocess.run(cmd, shell=True, check=True)
    # import the sam
    valid_candidates = []
    proportion_reference_covered = {}
    with pysam.AlignmentFile(output_file, "r") as sam_file:
        for read in sam_file.fetch():
            if read.is_unmapped:
                continue
            matching_bases = 0
            total_length = sam_file.get_reference_length(read.reference_name)
            for op, length in read.cigartuples:
                if op == 7:
                    matching_bases += length
            proportion_matching = (matching_bases / total_length) * 100
            if proportion_matching >= 98:
                if read.reference_name not in proportion_reference_covered:
                    proportion_reference_covered[read.reference_name] = set()
                proportion_reference_covered[read.reference_name].add(read.query_name)
    for ref in proportion_reference_covered:
        if len(proportion_reference_covered[ref]) >= 5:
            valid_candidates.append(ref)
    return valid_candidates

def process_assemblies(assembly_file, genes_file, output_json, nanopore):
    if not os.path.exists(os.path.dirname(output_json)):
        os.mkdir(os.path.dirname(output_json))
    if assembly_file.endswith('.fasta') or assembly_file.endswith('.fna') or assembly_file.endswith('.fa'):
        output_sam = os.path.join(os.path.dirname(output_json), f"{os.path.basename(assembly_file).replace('.fasta', '').replace('.fna', '').replace('.fa', '')}.sam")
        print(f"Processing {assembly_file}...")
        run_minimap2(assembly_file, genes_file, output_sam)
        hits = parse_sam(output_sam)
        filtered_hits = filter_longest_hits(hits)
        if nanopore is not None:
            additional_genes = supplement_with_genes_from_reads(genes_file, nanopore, output_sam.replace(".sam", ".supplement.sam"))
        else:
            additional_genes = []
        if ".json" in output_json:
            write_json(filtered_hits, output_json, assembly_file, additional_genes)
            write_gff(filtered_hits, output_json.replace(".json", ".gff"), assembly_file)
        if ".gff" in output_json:
            write_gff(filtered_hits, output_json.replace(".json", ".gff"), assembly_file)
        print(f"Finished processing {assembly_file}. Results saved to {output_json}.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Map genes to assemblies and generate GFF files of alleles with 95% similarity.")
    parser.add_argument("assembly_path", help="Assembly in FASTA format.")
    parser.add_argument("genes_file", help="FASTA file of genes to map to the assemblies.")
    parser.add_argument("output_file", help="Output file path for JSON.")
    parser.add_argument(
            "--nanopore",
            dest="nanopore",
            required=None,
    )
    args = parser.parse_args()
    process_assemblies(args.assembly_path, args.genes_file, args.output_file, args.nanopore)

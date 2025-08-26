import pandas as pd
import os

gene_list = "/hps/nobackup/iqbal/dander/pandora_systematic_evaluation/pandora_assessment_with_AMR/0.8_0/summary_stats.k15.w14_combined_plot.fp_genes.txt"
panaroo_pa = "/hps/nobackup/iqbal/dander/pandora_systematic_evaluation/pandora_assessment_with_AMR/0.8_0/panaroo_output/gene_presence_absence.csv"

with open(gene_list) as i:
    genes = set(i.read().split("\n"))

pa = pd.read_csv(panaroo_pa)
annotations = ["Gene name\tAnnotation"]
for index, row in pa.iterrows():
    if row["Gene"] in genes:
        annotations.append(f"{row['Gene']}\t{row['Annotation']}")
with open(gene_list.replace(".fp_genes.txt", ".fp_genes_annotations.txt"), "w") as o:
    o.write("\n".join(annotations))
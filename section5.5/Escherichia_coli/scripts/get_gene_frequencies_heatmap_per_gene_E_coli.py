import glob
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import seaborn as sns
import sys

sys.setrecursionlimit(50000)

amira_outputs = glob.glob("evaluation_results/Amira_output.v0.9.3/*/amira_results.tsv")

def apply_rules(gene):
    gene = gene.replace("'", "")
    if "blaCTX-M" in gene:
        gene = "blaCTX-M"
    if "blaNDM" in gene:
        gene = "blaNDM"
    if "blaOXA" in gene:
        gene = "blaOXA"
    if "aac(6)-Ib" in gene:
        gene = "aac(6)-Ib"
    if "blaEC" in gene:
        gene = "blaEC"
    if "oqx" in gene or "Oqx" in gene:
        gene = "oqx"
    if "blaTEM" in gene:
        gene = "blaTEM"
    if "blaOXY" in gene:
        gene = "blaOXY"
    if "fosA" in gene:
        gene = "fosA"
    if "blaCMY" in gene:
        gene = "blaCMY"
    if "aadA" in gene:
        gene = "aadA"
    if "arr-" in gene:
        gene = "arr-"
    if "dfrA" in gene:
        gene = "dfrA"
    if "rmtB" in gene:
        gene = "rmtB"
    if "aac(3)-VI" in gene:
        gene = "aac(3)-VI"
    if "aac(3)-II" in gene and "aac(3)-III" not in gene:
        gene = "aac(3)-II"
    if "aac(3)-III" in gene:
        gene = "aac(3)-III"
    if "blaSHV" in gene:
        gene = "blaSHV"
    if "qnrS" in gene:
        gene = "qnrS"
    if "aac(3)-I" in gene and "aac(3)-II" not in gene and "aac(3)-III" not in gene:
        gene = "aac(3)-I"
    if "blaKPC" in gene:
        gene = "blaKPC"
    if "mcr-5" in gene:
        gene = "mcr-5"
    if "mcr-1" in gene:
        gene = "mcr-1"
    if "mcr-2" in gene:
        gene = "mcr-2"
    if "qnrA" in gene:
        gene = "qnrA"
    if "qnrB" in gene:
        gene = "qnrB"
    if "cmlA" in gene:
        gene = "cmlA"
    if "aph(3)-I" in gene and "aph(3)-II" not in gene and "aph(3)-III" not in gene:
        gene = "aph(3)-I"
    if "aph(3)-II" in gene and "aph(3)-III" not in gene:
        gene = "aph(3)-II"
    if "aph(3)-III" in gene:
        gene = "aph(3)-III"
    if "aph(3)-VI" in gene:
        gene = "aph(3)-VI"
    if "aph(4)-I" in gene and "aph(4)-II" not in gene and "aph(4)-III" not in gene:
        gene = "aph(4)-I"
    if "aph(4)-II" in gene and "aph(4)-III" not in gene:
        gene = "aph(4)-II"
    if "aph(4)-III" in gene:
        gene = "aph(4)-III"
    if "aph(6)-I" in gene and "aph(6)-II" not in gene and "aph(6)-III" not in gene:
        gene = "aph(6)-I"
    if "aph(6)-II" in gene and "aph(6)-III" not in gene:
        gene = "aph(6)-II"
    if "aph(6)-III" in gene:
        gene = "aph(6)-III"
    if "qacE" in gene:
        gene = "qacE"
    if "blaLAP" in gene:
        gene = "blaLAP"
    if "aac(6)-I" in gene and "aac(6)-II" not in gene and "aac(6)-III" not in gene:
        gene = "aac(6)-I"
    if "aac(6)-II" in gene and "aac(6)-III" not in gene:
        gene = "aac(6)-II"
    if "aac(6)-III" in gene:
        gene = "aac(6)-III"
    if "blaDHA" in gene:
        gene = "blaDHA"
    if "qepA" in gene:
        gene = "qepA"
    if "blaIMI" in gene:
        gene = "blaIMI"
    if "ant(2)-I" in gene:
        gene = "ant(2)-I"
    if "ant(3)-I" in gene:
        gene = "ant(3)-I"
    if "ant(9)-I" in gene:
        gene = "ant(9)-I"
    if "blaACC" in gene:
        gene = "blaACC"
    if "blaACT" in gene:
        gene = "blaACT"
    if "blaCARB" in gene:
        gene = "blaCARB"
    if "blaGES" in gene:
        gene = "blaGES"
    if "blaIMP" in gene:
        gene = "blaIMP"
    if "blaGES" in gene:
        gene = "blaGES"
    if "blaVIM" in gene:
        gene = "blaVIM"
    if "toprJ" in gene:
        gene = "toprJ"
    if "blaVEB" in gene:
        gene = "blaVEB"
    if "blaMIR" in gene:
        gene = "blaMIR"
    if "arr" in gene:
        gene = "arr"
    if "mphK" in gene:
        gene = "mph(K)"
    if "satA" in gene:
        gene = "satA"
    if "blaOKP" in gene:
        gene = "blaOKP"
    if "blaHER" in gene:
        gene = "blaHER"
    if "qacG" in gene:
        gene = "qacG"
    if "blaMUN" in gene:
        gene = "blaMUN"
    if "blaMAL" in gene:
        gene = "blaMAL"
    if "blaORN" in gene:
        gene = "blaORN"
    if "blaPME" in gene:
        gene = "blaPME"
    if "catB" in gene:
        gene = "catB"
    if "cfxA" in gene:
        gene = "cfxA"
    if "tmexD" in gene:
        gene = "tmexD"
    if "blaPDC" in gene:
        gene = "blaPDC"
    if "mcr-3" in gene:
        gene = "mcr-3"
    if "mcr-8" in gene:
        gene = "mcr-8"
    if "catA" in gene:
        gene = "catA"
    if "qnrD" in gene:
        gene = "qnrD"
    if "blaLEN" in gene:
        gene = "blaLEN"
    if "blaPDC" in gene:
        gene = "blaPDC"
    if "blaPDC" in gene:
        gene = "blaPDC"
    if "blaSCO" in gene:
        gene = "blaSCO"
    if "rmtF" in gene:
        gene = "rmtF"
    if "blaFOX" in gene:
        gene = "blaFOX"
    if "blaADC" in gene:
        gene = "blaADC"
    if "blaFRI" in gene:
        gene = "blaFRI"
    if "blaCMH" in gene:
        gene = "blaCMH"
    if "blaSFO" in gene:
        gene = "blaSFO"
    if "cfiA" in gene:
        gene = "cfiA"
    if "cepA" in gene:
        gene = "cepA"
    if "blaGMA" in gene:
        gene = "blaGMA"
    if "ant(6)-I" in gene and not "ant(6)-II" in gene:
        gene = "ant(6)-I"
    return gene

# load all the genes we are interested in
with open("AMR_alleles_unified.fa") as i:
    allele_rows = i.read().split(">")[1:]
reference_genes = set()
ariba_mapping = {}
for r in allele_rows:
    if r != "":
        amira_allele, reference_allele = r.split("\n")[0].split(";")
        reference_genes.add(apply_rules(reference_allele.split(".NG")[0]))

total_counts = {}
count_data = {"sample": [], "gene": [], "value": []}
for a in tqdm(amira_outputs):
    sample = os.path.basename(os.path.dirname(a))
    amrfinder = os.path.join("evaluation_results/AMR_finder_plus_results.flye_v2.9.3_nanopore_only_assemblies", sample, "AMR_finder_plus_results.gff")
    if not os.path.exists(amrfinder):
        continue
    amira_content = pd.read_csv(a, sep="\t")
    amrfinder_content = pd.read_csv(amrfinder, sep="\t")
    # process amira results
    unique_amira_genes = {}
    for index, row in amira_content.iterrows():
        gene_name = apply_rules(row["Determinant name"])
        if gene_name not in unique_amira_genes:
            unique_amira_genes[gene_name] = 0
        unique_amira_genes[gene_name] += 1
    # process amrfinder results
    unique_amrfp_genes = {}
    for index, row in amrfinder_content.iterrows():
        if row["% Coverage of reference sequence"] < 85:
            continue
        gene_name = row["Gene symbol"]
        if row["Element subtype"] == "POINT":
            continue
        gene_name = apply_rules(gene_name)
        if gene_name not in reference_genes:
            continue
        if gene_name not in unique_amrfp_genes:
            unique_amrfp_genes[gene_name] = 0
        unique_amrfp_genes[gene_name] += 1
    # detect presence or absence
    for g in set(unique_amira_genes.keys()).union(unique_amrfp_genes):
        amrfp_count = unique_amrfp_genes.get(g, 0)
        amira_count = unique_amira_genes.get(g, 0)
        count_data["sample"].append(sample)
        count_data["gene"].append(g)
        count_data["value"].append(amira_count - amrfp_count)

# Pivot data for heatmap visualization
data_df = pd.DataFrame(count_data).pivot(index="sample", columns="gene", values="value").fillna(0)
print(data_df)
data_df.to_csv("evaluation_results/gene_count_data.tsv", sep="\t")
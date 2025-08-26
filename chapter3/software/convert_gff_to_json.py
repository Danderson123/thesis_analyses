import glob
import json
import os

gff_dir = "/hps/nobackup/iqbal/dander/thesis_figures/correction_assessment/assessment_results/truth_gffs"
files = glob.glob(os.path.join(gff_dir, "*.gff"))
outdir = "/hps/nobackup/iqbal/dander/thesis_figures/correction_assessment/assessment_results/truth_gff_jsons"

def parse_gff(gff_file):
    with open(gff_file) as i:
        gff_content, _ = i.read().split("##FASTA\n")
    gff_lines = {}
    for row in gff_content.split("\n"):
        if "Name=" in row:
            contig = row.split("\t")[0]
            if contig not in gff_lines:
                gff_lines[contig] = []
            gff_lines[contig].append(row.split("\t")[6] + row.split("Name=")[1])
    return gff_lines

if not os.path.exists(outdir):
    os.mkdir(outdir)
for f in files:
    gff_lines = parse_gff(f)
    with open(os.path.join(outdir, os.path.basename(f).replace(".gff", ".json")), "w") as o:
        o.write(json.dumps(gff_lines))

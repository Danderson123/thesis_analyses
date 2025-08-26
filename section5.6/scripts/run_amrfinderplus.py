import subprocess
import glob
import os

metaflye = glob.glob("metaflye_assemblies/*/assembly.fasta")
flye_out = "AMRFinderPlus_metaflye"
for a in metaflye:
    out = os.path.join(flye_out, os.path.basename(os.path.dirname(a)))
    if not os.path.exists(out):
        command = f"mkdir -p {flye_out} && amrfinder -d 2024-07-22.1 -n {a} -o {out} -t 4 --plus --organism Escherichia"
        subprocess.run(command, shell=True, check=True)

metaflye = glob.glob("nanomdbg_assemblies/*/contigs.fasta.gz")
mdgb_out = "AMRFinderPlus_nanomdbg"
for a in metaflye:
    out = os.path.join(mdgb_out, os.path.basename(os.path.dirname(a)))
    if not os.path.exists(out):
        command = f"mkdir -p {mdgb_out} && amrfinder -d 2024-07-22.1 -n {a} -o {out} -t 4 --plus --organism Escherichia"
        subprocess.run(command, shell=True, check=True)

flye = glob.glob("flye_assemblies/*/assembly.fasta")
flye_out = "AMRFinderPlus_flye"
for a in flye:
    out = os.path.join(flye_out, os.path.basename(os.path.dirname(a)))
    if not os.path.exists(out):
        command = f"mkdir -p {flye_out} && amrfinder -d 2024-07-22.1 -n {a} -o {out} -t 4 --plus --organism Escherichia"
        subprocess.run(command, shell=True, check=True)
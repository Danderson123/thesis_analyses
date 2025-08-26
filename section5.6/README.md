
# Overview

The scripts used in Section 5.6 of my thesis.

# Installation

The pipeline assumes the Amira singularity container is available in the directory and is called `amira.v0.11.0.img`. This can be installed [here](https://github.com/Danderson123/amira). You also need to have conda and singularity installed. All other dependencies can be installed via conda with:

```{bash}
conda env create -f mixed_strain_evaluation.yaml && conda activate mixed_strain_evaluation
```

The pipeline assumes the panRG is available in this directory and it can be downloaded from [here](https://figshare.com/ndownloader/files/56377289).

The pipeline assumes you have downloaded the truth evaluation data and it is available in the directory above this one in a subdirectory called `amira_paper_data`. The reference assemblies can be downloaded from [here](https://figshare.com/ndownloader/articles/28958849/versions/1). You will also need to download the Illumina and Nanopore reads and place them in subdirectories called `amira_paper_data/illumina_reads` and `amira_paper_data/nanopore_reads` respectively.

You will also need to download the AMRFinderPlus database (`v2024-07-22.1`) and have it available in this directory.

# Running

The reads can be downsampled and pooled using:
```{bash}
python3 scripts/downsample_reads.py && python3 scripts/merge_fastqs.py
```
The Flye assemblies can be generated with:
```{bash}
python3 scripts/flye_assemble.py
```
The metaFlye assemblies can be generated with:
```{bash}
python3 scripts/metaflye_assemble.py
```
The nanoMDBG assemblies can be generated with:
```{bash}
python3 scripts/nanomdbg_assemble.py
```
AMRFinderPlus can then be run on the assemblies with:
```{bash}
python3 scripts/run_amrfinderplus.py
```
Amira can be run with:
```{bash}
python3 scripts/run_amira.py
```

# Generating results

The final plot can be generated with this command after running the simulations for all read lengths:
```{bash}
python3 scripts/plot_recall_and_precision.py
```
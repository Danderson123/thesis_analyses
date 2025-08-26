
# Overview

The scripts used in Section 5.7 of my thesis.

# Installation

The pipeline assumes the Amira singularity container is available in the directory and is called `amira.v0.11.0.img`. This can be installed [here](https://github.com/Danderson123/amira). You also need to have conda and singularity installed. All other dependencies can be installed via conda with:

```{bash}
conda env create -f metagenome_mock_evaluation.yaml && conda activate metagenome_mock_evaluation
```

The pipeline assumes the ESKAPE pathogen panRG is available in this directory and it can be downloaded from [here](https://figshare.com/ndownloader/files/56861591).

The pipeline assumes you have downloaded the ESKAPE mock data from Purushothaman *et al.*(2024) and it is available in this directory in a subdirectory called `ESKAPE_mock_data`. It is available from [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB75510).

You will also need to download the AMRFinderPlus database (`v2024-07-22.1`) and have it available in this directory.

# Running

The metaFlye assemblies can be generated and AMRFinderPlus run on the assemblies with:
```{bash}
python3 scripts/metaflye_assemble.py
```
The nanoMDBG assemblies can be generated and AMRFinderPlus run on the assemblies with:
```{bash}
python3 scripts/nanomdbg_assemble.py
```
Amira can be run with:
```{bash}
python3 scripts/run_amira.py
```

# Generating results

The final plot (of just the first replicate for each DNA extraction kit) can be generated with this command after running the simulations for all read lengths:
```{bash}
python3 scripts/plot_recall_and_precision.py
```

Assuming that a large number of long-read FASTQs are available in a subdirectory called `Escherichia_coli/nanopore_reads`, the linear regression of read length SD against mean length can be run with:
```{bash}
python3 scripts/estimate_badread_params.py && python3 scripts/plot_mean_against_sd_read_lengths.py
```

After running all of the simulations, the simulated read lengths can be plotted with:
```{bash}
python3 scripts/plot_read_length_distributions.py
```
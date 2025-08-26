# Chapter 2 evaluation

The pipelines and scripts used in Chapter 2 of my thesis

## Installation

Install the dependencies with conda:

```{bash}
conda env create -f pandora_evaluation.yaml && conda activate pandora_evaluation
```

## Data

The data used in the pandora paper must be in the directory above this one in a directory called `pandora_data`. It can be downloaded from [here](https://ftp.ebi.ac.uk/pub/software/pandora/2021/data_no_fast5s.tar).

The data used in the amira paper must be in the directory above this one (including the reference assemblies and the nanopore reads) in a directory called `amira_paper_data`. Instructions to gather this data are available from [here](https://www.biorxiv.org/content/10.1101/2025.05.16.654303v2).

The bakta database `v5.0` must be in the directory above this and called `bakta_db`.

## Running

To recreate figures 2.1 and 2.2, run:
```{bash}
snakemake --cores <CORES> --use-conda --conda-frontend conda --snakefile Snakefile --nolock
```

To recreate figure 2.3, run:
```{bash}
snakemake --cores <CORES> --use-conda --conda-frontend conda --snakefile Snakefile_filtered --nolock
```

To recreate figure 2.4, run:
```{bash}
snakemake --cores <CORES> --use-conda --conda-frontend conda --snakefile Snakefile_filtered_with_test --nolock
```
and then run:
```{bash}
python3 software/plot_combined_test_sample_results.py
```

The final panRG-construction pipeline (figure 2.5) is available for download from [here](https://github.com/Danderson123/Amira_panRG_pipeline).
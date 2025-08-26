# Chapter 3 evaluation

The pipelines and scripts used in Chapter 3 of my thesis

## Installation

The pipeline assumes the Amira singularity container is available in the directory and is called `amira.v0.11.0.img`. This can be installed [here](https://github.com/Danderson123/amira). You also need to have conda and singularity installed. All other dependencies can be installed via conda with:

```{bash}
conda env create -f amira_correction_evaluation.yaml && conda activate amira_correction_evaluation
```

## Data

The data used in the pandora paper must be in the directory above this one in a directory called `pandora_data`. It can be downloaded from [here](https://ftp.ebi.ac.uk/pub/software/pandora/2021/data_no_fast5s.tar).

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
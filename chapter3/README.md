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

To run the pipeline:
```{bash}
snakemake --cores <CORES> --use-conda --conda-frontend conda --snakefile Snakefile --nolock
```

To recreate figure 3.3, run:
```{bash}
python3 software/plot_correction_trajectories.py
```
To recreate figure 3.5, run:
```{bash}
python3 software/align_AMR_reads_thesis_figure.py
```
To recreate figures 3.6, 3.7 and A.1 run:
```{bash}
python3 software/plot_true_context_unitig_alignments.py
```
To recreate figure 3.8 run:
```{bash}
python3 software/align_AMR_reads.py
```

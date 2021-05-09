Snakemake pipeline for reproducible and automated analysis of RNA-seq datasets deposited on GEO
================
Erik Zhivkoplias

May 08, 2021



## Introduction

TBD

## Snakemake pipeline execution

Snakemake is a workflow management system that ensures reproducible and scalable data analysis, highly populary in bioinformatics. It requires Python 3 and can be easily installed via Anaconda cloud service

### Download and install Conda

Download and install the latest Anaconda / Miniconda for your operation system. Example for Linux is provided below:

``` bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source .bashrc
```

### Clone repo and install conda environment

First clone this github repo:

``` bash
git clone https://github.com/zhivkoplias/work-in-progress/RNA_seq_RSEM_pipeline
cd RNA_seq_RSEM_pipeline
```

Now proceed with conda environment installation:

``` bash
conda create --name snakemake_env --file snakemake_env.yml
conda activate snakemake_env

```

### Run the pipeline :)

``` bash
snakemake -s Snakefile_pre --use-conda --conda-prefix test_salmon --cores 16; rm -r .snakemake/; snakemake -s Snakefile_GEO --use-conda --conda-prefix test_salmon --cores 16

```



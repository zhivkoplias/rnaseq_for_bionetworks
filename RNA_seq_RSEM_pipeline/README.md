# Snakemake pipeline for reproducible and automated analysis of RNA-seq datasets using RSEM

This repository contains a Snakemake workflow for processing RNA-seq data, starting from SRA accession numbers and producing gene expression count matrices. The pipeline is designed to be reproducible and scalable, handling data download, quality control, trimming, alignment, and quantification automatically.

The workflow executes the following sequence of tasks for each sample:
1.  **Data Retrieval**: Downloads single-end FASTQ data from the NCBI Sequence Read Archive (SRA).
2.  **Quality Control**: Runs FastQC to assess the quality of the raw reads.
3.  **Read Trimming**: Uses Trimmomatic to remove low-quality bases and poly-A tails.
4.  **Alignment**: Aligns trimmed reads to a reference genome using Bowtie 2.
5.  **Expression Quantification**: Calculates gene and isoform expression levels from the alignments using RSEM.
6.  **Matrix Generation**: Merges the expression counts from all samples into final matrices suitable for downstream analysis.

---
## Snakemake pipeline execution

Snakemake is a workflow management system that ensures reproducible and scalable data analysis, highly popular in bioinformatics. It requires Python 3 and can be easily installed via the Conda package manager.

### 1. Download and install Conda

Download and install the latest Anaconda or Miniconda for your operating system. The example for Linux is provided below:

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### 2. Clone repo and set up environment

First, clone this GitHub repo:

```bash
git clone https://github.com/zhivkoplias/work-in-progress/RNA_seq_RSEM_pipeline
cd RNA_seq_RSEM_pipeline
```
This pipeline uses `--use-conda` to automatically manage all required software dependencies. You only need a base environment with Snakemake installed.

```bash
# Create a base environment for snakemake
conda create -n snakemake "snakemake>=5.3.0"
conda activate snakemake
```

### 3. Configure the pipeline

Before running, you must configure the paths and parameters.

**A. `config.yml`**:
All pipeline settings are controlled via `config.yml`. You must specify input files and can adjust tool parameters.

```yaml
# Required Paths
sample_file: "samples.tsv"
GSM_file: "GSM_info.tsv"

# Resource Settings
n_cores: 8

# Trimmomatic & Bowtie 2 Params...
# (see file for all options)
```

**B. `samples.tsv`**:
Create a `samples.tsv` file (or the name specified in `config.yml`). This file maps your desired sample names to their SRA accession numbers. It must contain two space-separated columns with no header:

**Example `samples.tsv`**:
```
SRR8941013 wildtype_rep1
SRR8941014 wildtype_rep2
SRR8941015 gadX_knockout_rep1
SRR8941016 gadX_knockout_rep2
```

**C. Reference Genome**:
This pipeline requires pre-built **Bowtie 2** and **RSEM** indices. Place these files in the `data/reference/` directory as specified in the `Snakefile`.

### 4. Run the pipeline :)

Execute the workflow from the repository's root directory. The `--use-conda` flag will automatically download and set up all required bioinformatics tools.

To run the full pipeline using the main `Snakefile`:
```bash
snakemake --use-conda --cores <number_of_cores>
```
*   `--use-conda`: Instructs Snakemake to manage software environments.
*   `--cores`: Specifies the total number of CPU cores to use.

If your workflow is split into two parts (e.g., fetching data then processing), you can use your original command structure:
```bash
# Example for a two-stage execution
snakemake -s Snakefile_pre --use-conda --cores 16; \
rm -rf .snakemake/; \
snakemake -s Snakefile_GEO --use-conda --cores 16
```

---
## Output Directory Structure

The pipeline will generate the following directory structure:

*   `data/processed/`: Contains downloaded (`.fastq.gz`), trimmed (`.trimmed.fastq.gz`), and aligned (`.bam`) files.
*   `results/counts/`: Contains RSEM expression results (`.genes.results`, `.isoforms.results`) for each sample.
*   `results/reports/fastqc/`: Contains FastQC quality reports for each raw FASTQ file.
*   `results/matrices/`: Contains the final merged expression matrices.
*   `logs/`: Contains log files for each step of the pipeline.


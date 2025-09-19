# RNASeq Analysis Pipelines for network inference

This repository is a collection of pipelines and scripts for processing and analyzing gene expression data so it can be used with tools and algorithms for network inference. It includes workflows for shRNA‑seq, standard RNA‑seq, and L1000 data. Each pipeline is self‑contained in its own directory.

---

## Pipelines in This Repository

Below is a summary of the analysis workflows available here. For detailed instructions on setup and execution, please refer to the `README.md` file located within each pipeline's respective folder.

### 1. shRNA-seq Perturbation Analysis Pipeline (`RNA_seq_Salmon_pipeline/`)

This pipeline analyzes shRNA-seq perturbation experiments from the ENCODE project. It uses a pseudo-alignment approach with Salmon to quantify transcript abundance and calculate log2 fold-change values.

**Overview:**
The workflow is divided into three main stages:
1.  **Data Acquisition**: Fetching metadata from the ENCODE portal and filtering for relevant experiment files (`.bam`).
2.  **Transcript Quantification**: Running Salmon on the filtered `.bam` files to get transcript-level counts (TPM).
3.  **Post-processing & FC Calculation**: Merging counts by gene, associating perturbations with controls, and calculating log2FC values to generate the final output matrices suitable for downstream modeling.

> For detailed instructions, please see the README inside the [`RNA_seq_Salmon_pipeline/`](./RNA_seq_Salmon_pipeline/) folder.

---

### 2. Snakemake RNA-seq Pipeline (`RNA_seq_RSEM_pipeline/`)

This is a reproducible and automated Snakemake workflow for processing standard RNA-seq datasets. It starts with raw data from the Sequence Read Archive (SRA) and uses a traditional alignment-based approach with Bowtie 2 and RSEM.

**Overview:**
The workflow executes the following sequence of tasks for each sample:
1.  **Data Retrieval**: Downloads FASTQ data from the NCBI Sequence Read Archive (SRA).
2.  **Quality Control**: Runs FastQC to assess raw read quality.
3.  **Read Trimming**: Uses Trimmomatic to remove low-quality bases.
4.  **Alignment**: Aligns trimmed reads to a reference genome using Bowtie 2.
5.  **Expression Quantification**: Calculates gene and isoform expression levels using RSEM.
6.  **Matrix Generation**: Merges expression counts from all samples into final matrices.

> For detailed instructions, please see the README inside the [`RNA_seq_RSEM_pipeline/`](./RNA_seq_RSEM_pipeline/) folder.

---

### 3. GCTX Data Parsing Wrappers (`GCTX_counts_L1000_project/`)

This directory contains a set of wrapper scripts to parse and handle GCTX data from the L1000 project, primarily using the `cmapPy` library.

**Overview:**

TBD


# shRNA-seq Perturbation Analysis Pipeline using Salmon

This repository contains a set of scripts to analyze shRNA-seq perturbation experiments from the ENCODE project for the K562 and HepG2 cell lines. The pipeline processes raw alignment files, quantifies transcript abundance using Salmon, and calculates log2 fold-change (log2FC) values to produce matrices suitable for downstream modeling.

## Overview

The workflow is divided into three main stages:
1.  **Data Acquisition**: Fetching metadata from the ENCODE portal and filtering for relevant experiment files (`.bam`).
2.  **Transcript Quantification**: Running Salmon on the filtered `.bam` files to get transcript-level counts (TPM).
3.  **Post-processing & FC Calculation**: Merging counts by gene, associating perturbations with controls, and calculating log2FC values to generate the final output matrices.

---

## Prerequisites

Before running the pipeline, ensure you have the following software installed:

*   **Python 3**: with `requests` library (`pip install requests`)
*   **Salmon**: For transcript quantification.
*   **BEDTools**: Specifically, the `bamToFastq` utility is required.
*   **R**: with the following libraries installed:
    *   `tximport`
    *   `dplyr`
    *   `stringr`
    *   `GenomicFeatures`
    *   `mygene`
*   **A reference transcriptome index for Salmon**: The pipeline was run using an index built from **GRCh38** (Gencode release 37). You must build this yourself using the reference transcriptome.

---

## Pipeline Steps

### 1. Data Acquisition

This step identifies and filters the required experiment and control files from the ENCODE database.

**Action**:
First, download two file lists from the ENCODE project portal:
1.  **Perturbation Experiments (shRNA RNA-seq, human)**
    > `https://www.encodeproject.org/search/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&perturbed=true&biosample_ontology.term_name=K562&assay_title=shRNA+RNA-seq&limit=all`
2.  **Control Experiments (shRNA RNA-seq, human, not perturbed)**
    > `https://www.encodeproject.org/search/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=K562&assay_title=shRNA+RNA-seq&perturbed=false&limit=all`

From the downloaded metadata, create two simple text files containing the download URLs of the `.bam` files: one for perturbations (`perturbations_raw.txt`) and one for controls (`controls_raw.txt`).

Next, run the `parseENCODE.py` script to filter these files based on specific criteria (assembly: GRCh38, output type: transcriptome alignments). This script also generates a mapping file between file IDs and experiment IDs.

**Usage**:
```bash
# Filter perturbation files
python parseENCODE.py \
    -f1 perturbations_raw.txt \
    -f2 ENCODE_perturbations_filtered.txt \
    -f3 ENCODE_perturbations_files2experiments.csv

# Filter control files
python parseENCODE.py \
    -f1 controls_raw.txt \
    -f2 ENCODE_controls_filtered.txt \
    -f3 ENCODE_controls_files2experiments.csv
```

Finally, use the filtered lists (`ENCODE_perturbations_filtered.txt`, `ENCODE_controls_filtered.txt`) to download the required `.bam` files to a local directory. **Note**: The total data size is approximately 2.5 TB.

### 2. Transcript Quantification

This step runs Salmon to quantify transcript abundance from the downloaded `.bam` files.

**Action**:
The `run_salmon_experiments.sh` script is designed for a SLURM-based HPC cluster. It iterates through each `.bam` file, converts it to FASTQ using `bamToFastq`, and then runs `salmon quant`.

**Key Script Logic**:
```bash
# For each BAM file...
BAM_FILE=$i
FASTQ_PREFIX="${BAM_FILE}.fq"

# 1. Convert BAM to paired-end FASTQ
bamToFastq -i $BAM_FILE -fq ${FASTQ_PREFIX}1.fq -fq2 ${FASTQ_PREFIX}2.fq

# 2. Run Salmon
salmon quant -i <path_to_human_index> -l A \
    -1 ${FASTQ_PREFIX}1.fq \
    -2 ${FASTQ_PREFIX}2.fq \
    --validateMappings \
    -o <output_directory>/${BAM_FILE_BASENAME}

# 3. Clean up FASTQ files
rm ${FASTQ_PREFIX}1.fq ${FASTQ_PREFIX}2.fq
```
**Usage**:
You will need to modify `run_salmon_experiments.sh` to point to your file locations and specify your Salmon index. This is a computationally intensive step that will take several days to complete.

### 3. Post-processing and Fold-Change Calculation

This step processes the raw Salmon output (`quant.sf` files) into final log2FC matrices.

**Action**:
The `run_salmon_counts_to_FC.Rmd` script performs the following key operations:
1.  **Loads Data**: Reads the `files2experiments` mapping tables and experiment metadata.
2.  **Imports Counts**: Uses `tximport` to load Salmon TPM values and aggregates them from the transcript level to the gene level.
3.  **Calculates log2FC**: For each perturbation experiment, it identifies the corresponding control experiment, calculates the log2 fold-change `log2((TPM_perturbation + pseudo_count) / (TPM_control + pseudo_count))`, and handles replicates by averaging.
4.  **Aggregates by Gene**: If multiple shRNA experiments target the same gene, their results are aggregated to create a final gene-by-gene matrix.
5.  **Generates Matrices**: Produces the final output files.

**Usage**:
Open `run_salmon_counts_to_FC.Rmd` in RStudio or an R environment and execute the code chunks. Ensure all file paths at the beginning of the script are updated to match your directory structure.

---

## Final Output

The pipeline generates two primary files:

1.  `ymatrix.csv`: The **data matrix (Y)**.
    *   **Rows**: Genes
    *   **Columns**: Perturbation experiments (shRNAs targeting a gene)
    *   **Values**: `log2FC` of gene expression resulting from the perturbation.

2.  `pmatrix.csv`: The **perturbation matrix (P)**.
    *   This is a design matrix where non-zero entries (0 or -1) indicate which gene is being directly targeted in a given perturbation experiment (column). It is structured to match the `Y` matrix.


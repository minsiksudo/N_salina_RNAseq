#!/bin/bash

# Define directories
RAW_READS_DIR="Nannochloropsis_20eax4_RNAseq"
QC_DIR="QC_results"

# Create output directory if it doesn't exist
mkdir -p $QC_DIR

# Activate the FastQC environment
conda activate fastqc_env

# Run FastQC
fastqc -o $QC_DIR $RAW_READS_DIR/*.fq.gz

# Deactivate the FastQC environment
conda deactivate

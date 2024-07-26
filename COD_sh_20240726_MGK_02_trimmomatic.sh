#!/bin/bash

# Define directories
RAW_READS_DIR="Nannochloropsis_20eax4_RNAseq"
FILTERED_DIR="filtered_reads"
ADAPTERS_DIR="adapters"

# Set threads for parallel computing
THREADS=40

# Create output directory if it doesn't exist
mkdir -p $FILTERED_DIR

# Activate the Trimmomatic environment
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate trimmomatic_env

# Run Trimmomatic
for file in $RAW_READS_DIR/*_1.fq.gz; do
    BASENAME=$(basename $file _1.fq.gz)
    trimmomatic PE -threads $THREADS \
        ${RAW_READS_DIR}/${BASENAME}_1.fq.gz ${RAW_READS_DIR}/${BASENAME}_2.fq.gz \
        ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz ${FILTERED_DIR}/${BASENAME}_1_unpaired.fq.gz \
        ${FILTERED_DIR}/${BASENAME}_2_paired.fq.gz ${FILTERED_DIR}/${BASENAME}_2_unpaired.fq.gz \
        ILLUMINACLIP:${ADAPTERS_DIR}/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36

done

# Deactivate the Trimmomatic environment
conda deactivate

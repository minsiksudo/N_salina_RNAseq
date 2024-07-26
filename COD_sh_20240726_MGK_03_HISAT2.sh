#!/bin/bash

# Define directories
FILTERED_DIR="filtered_reads"
ALIGN_DIR="aligned_reads"
REF_GENOME="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/GCA_004565275.1_ASM456527v1_genomic.fna"
HISAT2_INDEX="hisat2_index/nannochloropsis_salina"

# Set threads 
THREADS=40

# Create output directory if it doesn't exist
mkdir -p $ALIGN_DIR

#Activate the HISAT2 env

source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate hisat2_env

# Build HISAT2 index
hisat2-build $REF_GENOME $HISAT2_INDEX

# Align reads
for file in $FILTERED_DIR/*_1_paired.fq.gz; do
    BASENAME=$(basename $file _1_paired.fq.gz)
    hisat2 -p $THREADS -x $HISAT2_INDEX -1 ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz -2 ${FILTERED_DIR}/${BASENAME}_2_paired.fq.gz -S ${ALIGN_DIR}/${BASENAME}_aligned.sam
done

conda deactivate

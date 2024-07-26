#!/bin/bash

# Define the number of threads
THREADS=40

# Define directories
ALIGNED_DIR="aligned_reads"
FILTERED_DIR="filtered_reads"

# Create output directory if it doesn't exist
mkdir -p $FILTERED_DIR

# Activate the Samtools environment
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate samtools_env

# Run Samtools to filter aligned reads
for file in $ALIGNED_DIR/*.sam; do
    BASENAME=$(basename $file .sam)
    samtools view -@ $THREADS -bS $file | samtools sort -@ $THREADS -o ${FILTERED_DIR}/${BASENAME}_sorted.bam
    samtools index ${FILTERED_DIR}/${BASENAME}_sorted.bam
    samtools view -@ $THREADS -b -q 20 ${FILTERED_DIR}/${BASENAME}_sorted.bam > ${FILTERED_DIR}/${BASENAME}_filtered.bam
    samtools index ${FILTERED_DIR}/${BASENAME}_filtered.bam
done

# Deactivate the Samtools environment
conda deactivate

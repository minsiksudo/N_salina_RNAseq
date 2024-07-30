#!/bin/bash

# Define directories
FILTERED_DIR="filtered_reads"
COUNTS_DIR="counts"

# Path to annotation file
ANNOTATION_FILE="/path/to/annotation.gtf"

# Define the number of threads
THREADS=40

# Create directory if it doesn't exist
mkdir -p $COUNTS_DIR

# Generate count tables for each BAM file
for file in $FILTERED_DIR/*.bam; do
    BASENAME=$(basename $file .bam)
    featureCounts -a $ANNOTATION_FILE -o ${COUNTS_DIR}/${BASENAME}_counts.txt -t exon -g gene_id -T $THREADS $file
done

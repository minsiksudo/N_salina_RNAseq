#!/bin/bash

# Define variables
BAM_DIR="aligned_reads"  # Directory containing SAM/BAM files
GTF_FILE="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/genomic.gtf"  # Path to GTF file
OUTPUT_DIR="counts"  # Directory to store featureCounts results
THREADS=40  # Number of threads to use

# Create output directories if they do not exist
mkdir -p $OUTPUT_DIR

# Step 1: Convert SAM files to BAM (if SAM exists)
for SAM in $BAM_DIR/*.sam; do
    BAM=${SAM/.sam/.bam}
    if [ ! -f "$BAM" ]; then
        echo "Converting $SAM to BAM format..."
        samtools view -S -b $SAM > $BAM
    fi
done

# Step 2: Sort BAM files
for BAM in $BAM_DIR/*.bam; do
    SORTED_BAM=${BAM/.bam/_sorted.bam}
    if [ ! -f "$SORTED_BAM" ]; then
        echo "Sorting $BAM..."
        samtools sort $BAM -o $SORTED_BAM
    fi
done

# Step 3: Index sorted BAM files
for SORTED_BAM in $BAM_DIR/*_sorted.bam; do
    if [ ! -f "${SORTED_BAM}.bai" ]; then
        echo "Indexing $SORTED_BAM..."
        samtools index $SORTED_BAM
    fi
done

# Step 4: Run featureCounts on the sorted BAM files
for SORTED_BAM in $BAM_DIR/*_sorted.bam; do
    BASENAME=$(basename $SORTED_BAM _sorted.bam)
    COUNTS_FILE="${OUTPUT_DIR}/${BASENAME}_counts.txt"
    if [ ! -f "$COUNTS_FILE" ]; then
        echo "Running featureCounts on $SORTED_BAM..."
        featureCounts -a $GTF_FILE -o $COUNTS_FILE -T $THREADS $SORTED_BAM
    fi
done

echo "All steps completed successfully!"

#!/bin/bash

# Define the number of threads and memory per thread
THREADS=40
MEMORY_PER_THREAD="8G"  # Allocate 8 GB of memory per thread (adjust as needed)

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
    
    # Skip processing if the file is empty or does not exist
    if [[ ! -s $file ]]; then
        echo "Warning: Skipping empty or non-existent file $file"
        continue
    fi
    
    echo "Processing $file..."
    
    # Convert SAM to BAM, sort, and index
    samtools view -@ $THREADS -bS $file | \
    samtools sort -@ $THREADS -m $MEMORY_PER_THREAD -o ${FILTERED_DIR}/${BASENAME}_sorted.bam

    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to convert and sort $file"
        continue
    fi

    samtools index ${FILTERED_DIR}/${BASENAME}_sorted.bam

    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to index sorted BAM for $file"
        continue
    fi

    samtools view -@ $THREADS -b -q 20 ${FILTERED_DIR}/${BASENAME}_sorted.bam > ${FILTERED_DIR}/${BASENAME}_filtered.bam

    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to filter BAM for $file"
        continue
    fi

    samtools index ${FILTERED_DIR}/${BASENAME}_filtered.bam

    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to index filtered BAM for $file"
        continue
    fi

    echo "Finished processing $file"
done

# Deactivate the Samtools environment
conda deactivate

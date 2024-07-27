#!/bin/bash

# Define the number of threads
THREADS=40

# Define directories
FILTERED_DIR="filtered_reads"
ANNOTATED_DIR="annotated_reads"

# Create output directory if it doesn't exist
mkdir -p $ANNOTATED_DIR

# Activate the eggNOG-mapper environment
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_mapper_env

# Run EggNOG-mapper
for file in $FILTERED_DIR/*.bam; do
    BASENAME=$(basename $file _filtered.bam)
    emapper.py -i $file -o ${ANNOTATED_DIR}/${BASENAME} --cpu $THREADS
done

# Deactivate the eggNOG-mapper environment
conda deactivate

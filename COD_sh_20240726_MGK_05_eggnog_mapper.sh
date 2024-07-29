#!/bin/bash

# Define the number of threads
THREADS=40

# Define directories
FILTERED_DIR="filtered_reads"
ANNOTATED_DIR="annotated_reads"

# Path to custom EggNOG database directory
EGGNOG_DB_DIR="/mnt/4T_samsung/Dropbox/Sequencing archive/eggnog_db"

# Create output directory if it doesn't exist
mkdir -p $ANNOTATED_DIR

# Activate the eggNOG-mapper environment
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_mapper_env

# Run EggNOG-mapper with the custom database path
for file in $FILTERED_DIR/*.bam; do
    BASENAME=$(basename $file _filtered.bam)
    emapper.py -i $file -o ${ANNOTATED_DIR}/${BASENAME} --cpu $THREADS --data_dir $EGGNOG_DB_DIR
done

# Deactivate the eggNOG-mapper environment
conda deactivate

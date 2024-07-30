#!/bin/bash

# Define the number of threads
THREADS=40

# Define directories
FILTERED_DIR="filtered_reads"
FASTA_DIR="fasta_reads"
ANNOTATED_DIR="annotated_reads"

# Path to custom EggNOG database directory
EGGNOG_DB_DIR="/mnt/4T_samsung/Dropbox/Sequencing_archive/eggnog_db"

# Create directories if they don't exist
mkdir -p $FASTA_DIR
mkdir -p $ANNOTATED_DIR

# Activate the eggNOG-mapper environment
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_mapper_env

# Convert BAM files to FASTA format
for file in $FILTERED_DIR/*.bam; do
    BASENAME=$(basename $file .bam)
    samtools fasta $file > ${FASTA_DIR}/${BASENAME}.fasta
done

# Run EggNOG-mapper with the custom database path
for file in $FASTA_DIR/*.fasta; do
    BASENAME=$(basename $file .fasta)
    emapper.py --cpu $THREADS --output ${ANNOTATED_DIR}/${BASENAME} --data_dir "$EGGNOG_DB_DIR" -i $file --override --translate -m diamond
done

# Deactivate the eggNOG-mapper environment
conda deactivate

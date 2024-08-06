#!/bin/bash

# Define directories and number of threads
FILTERED_DIR="filtered_reads"
FASTA_DIR="fasta_reads"
ANNOTATED_DIR="annotated_reads"
THREADS=40
EGGNOG_DB_DIR="/mnt/4T_samsung/Dropbox/Sequencing_archive/eggnog_db"

# Create directories if they don't exist
mkdir -p $FASTA_DIR
mkdir -p $ANNOTATED_DIR

# Activate the eggNOG-mapper environment
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_mapper_env

# Create a list of already annotated files (without extensions)
ls $ANNOTATED_DIR/*.emapper.annotations | sed 's/\.[^.]*$//' | xargs -n 1 basename > annotated_files.txt

# Convert BAM files to FASTA format
for file in $FILTERED_DIR/*.bam; do
    BASENAME=$(basename $file .bam)
    FASTA_FILE="${FASTA_DIR}/${BASENAME}.fasta"
    
    # Convert BAM to FASTA if not already in FASTA directory
    if [ ! -f $FASTA_FILE ]; then
        samtools fasta $file > $FASTA_FILE
    fi
done

# Run EggNOG-mapper with the custom database path
for file in $FASTA_DIR/*.fasta; do
    BASENAME=$(basename $file .fasta)

    # Skip files that have already been annotated
    if ! grep -q "^$BASENAME$" annotated_files.txt; then
        emapper.py --cpu $THREADS --output ${ANNOTATED_DIR}/${BASENAME} --data_dir "$EGGNOG_DB_DIR" -i $file --override --translate -m diamond
    fi
done

# Clean up intermediate FASTA files after annotation
echo "Cleaning up intermediate FASTA files..."
rm -f ${FASTA_DIR}/*.fasta

# Deactivate the eggNOG-mapper environment
conda deactivate

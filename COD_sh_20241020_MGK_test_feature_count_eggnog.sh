#!/bin/bash

# Set directories
FASTA_DIR="fasta_reads"
ANNOTATED_DIR="annotated_reads"
COUNTS_DIR="counts"

# Create output directories if they don't exist
mkdir -p $COUNTS_DIR

# Activate the environment if needed (this depends on whether featureCounts or eggnog_mapper requires a conda environment)
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate eggnog_mapper_env  # if necessary

# Define the sample to process (replace with your actual sample)
BASENAME="TN1506R0382"
ANNOTATION_FILE="$ANNOTATED_DIR/${BASENAME}_aligned_filtered.emapper.annotations"
FASTA_FILE="$FASTA_DIR/${BASENAME}.fasta"

# Check if the FASTA and annotation files exist
if [ ! -f "$FASTA_FILE" ]; then
    echo "FASTA file for $BASENAME not found. Exiting..."
    exit 1
fi

if [ ! -f "$ANNOTATION_FILE" ]; then
    echo "Annotation file for $BASENAME not found. Exiting..."
    exit 1
fi

# Function to convert EggNOG annotations to SAF format (in-memory using a process substitution)
# Example parsing: Adjust this based on the structure of your EggNOG annotations
# This assumes the EggNOG annotation file contains: GeneID, Chr, Start, End, Strand
convert_to_saf() {
    echo -e "GeneID\tChr\tStart\tEnd\tStrand"
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' "$ANNOTATION_FILE"
}

# Run featureCounts using process substitution to avoid creating intermediate files
featureCounts -F SAF -a <(convert_to_saf) -o "$COUNTS_DIR/${BASENAME}_counts.txt" "$FASTA_FILE"

if [ $? -ne 0 ]; then
    echo "Error running featureCounts for $BASENAME"
    exit 1
fi

echo "Feature counting completed for $BASENAME."

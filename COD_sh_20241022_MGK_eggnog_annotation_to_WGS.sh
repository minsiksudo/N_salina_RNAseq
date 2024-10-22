#!/bin/bash

# Script: run_eggnog_annotation.sh
# Purpose: Automate the process of running EggNOG-mapper and generating functional annotations for Nannochloropsis salina

# Directories and input/output files


# Define variables
PROTEIN_FILE="/mnt/4T_samsung/Dropbox/Sequencing_archive/2024_01_02_KNK_NSalina/ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/protein.faa"
OUTPUT_PREFIX="nannochloropsis_annotation"
EGGNOG_DB_DIR="/mnt/4T_samsung/Dropbox/Sequencing_archive/eggnog_db"
THREADS=40

# EggNOG-mapper command
echo "Running EggNOG-mapper on ${PROTEIN_FILE} with ${THREADS} threads..."

	source /home/bagel/miniconda3/etc/profile.d/conda.sh
	conda init
	conda activate eggnog_mapper_py3_env

emapper.py -i ${PROTEIN_FILE} \
           --output ${OUTPUT_PREFIX} \
           --cpu ${THREADS} \
           -m diamond \
           --data_dir ${EGGNOG_DB_DIR} \
           --dmnd_db ${EGGNOG_DB_DIR}/eggnog_proteins.dmnd \
           --go_evidence experimental \
           --tax_scope auto \
           --decorate_gff /mnt/4T_samsung/Dropbox/Sequencing_archive/2024_01_02_KNK_NSalina/ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/genomic_eggnog.gff \
	   --override \
           --output_dir ./eggnog_output

if [ $? -eq 0 ]; then
    echo "EggNOG-mapper completed successfully."
else
    echo "EggNOG-mapper encountered an error."
    exit 1
fi

conda deactivate

# Check the output files
echo "Annotation results saved in: ./eggnog_output/${OUTPUT_PREFIX}.emapper.annotations"

# TODO: Convert EggNOG annotations to GTF format
# This step needs specific logic to parse the EggNOG annotations and map them to GTF format.
# Depending on how you'd like the hierarchy information structured, custom scripts would be required here.
# For example, you can parse the .annotations file and assign functional descriptions to genes based on their positions.

# Placeholder message for GTF conversion
echo "Next step: Parse the EggNOG annotations and convert them to GTF format (this part is manual or needs a custom script)."

exit 0

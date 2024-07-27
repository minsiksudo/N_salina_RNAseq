#!/bin/bash

# Define directories
FILTERED_DIR="filtered_reads"
ALIGN_DIR="aligned_reads"
REF_GENOME="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/GCA_004565275.1_ASM456527v1_genomic.fna"
HISAT2_INDEX="hisat2_index/nannochloropsis_salina"

# Set threads 
THREADS=20

# Create output directory if it doesn't exist
mkdir -p $ALIGN_DIR

# Activate the HISAT2 environment
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate hisat2_env

# Build HISAT2 index
echo "Building HISAT2 index..." | tee build_index.log
if ! hisat2-build $REF_GENOME $HISAT2_INDEX &>> build_index.log; then
    echo "Error building HISAT2 index. Check build_index.log for details." | tee -a build_index.log
    conda deactivate
    exit 1
fi

# Align reads
for file in $FILTERED_DIR/*_1_paired.fq.gz; do
    BASENAME=$(basename $file _1_paired.fq.gz)
    echo "Aligning ${BASENAME}..." | tee -a align.log
    
    # Align with HISAT2
    if ! hisat2 -p $THREADS -x $HISAT2_INDEX -1 ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz -2 ${FILTERED_DIR}/${BASENAME}_2_paired.fq.gz -S ${ALIGN_DIR}/${BASENAME}_aligned.sam &> ${ALIGN_DIR}/${BASENAME}_align.log; then
        echo "Error aligning ${BASENAME}. Check ${ALIGN_DIR}/${BASENAME}_align.log for details." | tee -a align.log
    else
        echo "Finished aligning ${BASENAME}." | tee -a align.log
    fi
done

# Deactivate the HISAT2 environment
conda deactivate

echo "HISAT2 alignment process completed." | tee -a align.log

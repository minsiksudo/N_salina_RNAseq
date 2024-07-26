#!/bin/bash

# Define directories
RAW_READS_DIR="Nannochloropsis_20eax4_RNAseq"
FILTERED_DIR="filtered_reads"
ANNOTATED_DIR="annotated_reads"
SUMMARY_FILE="results/read_counts_summary.csv"

# Create results directory if it doesn't exist
mkdir -p results

# Initialize summary file
echo "sample,raw_reads,high_quality_reads,reference_genome_aligned_reads,annotated_reads,proportion_contaminants,proportion_annotated" > $SUMMARY_FILE

# Summarize the results
for file in $RAW_READS_DIR/*_1.fq.gz; do
    BASENAME=$(basename $file _1.fq.gz)
    RAW_READS=$(zcat ${RAW_READS_DIR}/${BASENAME}_1.fq.gz | wc -l)
    HQ_READS=$(zcat ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz | wc -l)
    ALIGN_READS=$(samtools view -c ${FILTERED_DIR}/${BASENAME}_filtered.bam)
    ANNOT_READS=$(grep -c ">" ${ANNOTATED_DIR}/${BASENAME}_annotation.emapper.annotations)
    PROP_CONTAM=$(echo "scale=2; 100 - ($ALIGN_READS / $RAW_READS * 100)" | bc)
    PROP_ANNOT=$(echo "scale=2; $ANNOT_READS / $RAW_READS * 100" | bc)
    echo "$BASENAME,$RAW_READS,$HQ_READS,$ALIGN_READS,$ANNOT_READS,$PROP_CONTAM,$PROP_ANNOT" >> $SUMMARY_FILE
done

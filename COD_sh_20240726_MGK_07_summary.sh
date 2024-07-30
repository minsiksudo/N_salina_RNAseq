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

    # Check if files exist
    if [ ! -f ${RAW_READS_DIR}/${BASENAME}_1.fq.gz ]; then
        echo "Raw reads file ${RAW_READS_DIR}/${BASENAME}_1.fq.gz not found."
        continue
    fi
    if [ ! -f ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz ]; then
        echo "Filtered reads file ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz not found."
        continue
    fi
    if [ ! -f ${FILTERED_DIR}/${BASENAME}_filtered.bam ]; then
        echo "Filtered BAM file ${FILTERED_DIR}/${BASENAME}_filtered.bam not found."
        continue
    fi
    if [ ! -f ${ANNOTATED_DIR}/${BASENAME}_annotation.emapper.annotations ]; then
        echo "Annotation file ${ANNOTATED_DIR}/${BASENAME}_annotation.emapper.annotations not found."
        continue
    fi

    RAW_READS=$(($(zcat ${RAW_READS_DIR}/${BASENAME}_1.fq.gz | wc -l) / 4))
    HQ_READS=$(($(zcat ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz | wc -l) / 4))
    ALIGN_READS=$(samtools view -c ${FILTERED_DIR}/${BASENAME}_filtered.bam)
    ANNOT_READS=$(grep -c ">" ${ANNOTATED_DIR}/${BASENAME}_annotation.emapper.annotations)
    
    if [ $RAW_READS -eq 0 ]; then
        PROP_CONTAM="N/A"
        PROP_ANNOT="N/A"
    else
        PROP_CONTAM=$(echo "scale=2; 100 - ($ALIGN_READS / $RAW_READS * 100)" | bc)
        PROP_ANNOT=$(echo "scale=2; $ANNOT_READS / $RAW_READS * 100" | bc)
    fi
    
    echo "$BASENAME,$RAW_READS,$HQ_READS,$ALIGN_READS,$ANNOT_READS,$PROP_CONTAM,$PROP_ANNOT" >> $SUMMARY_FILE
done

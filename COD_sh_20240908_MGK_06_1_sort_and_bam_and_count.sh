#!/bin/bash

# Define directories
RAW_READS_DIR="Nannochloropsis_20eax4_RNAseq"
FILTERED_DIR="filtered_reads"
ANNOTATED_DIR="annotated_reads"
SUMMARY_FILE="results/read_counts_summary.csv"
NUM_THREADS=40  # Define number of threads (adjust based on your system)

# Create results directory if it doesn't exist
mkdir -p results

# Initialize summary file
echo "sample,raw_reads,high_quality_reads,reference_genome_aligned_reads,annotated_reads,proportion_contaminants,proportion_annotated" > "$SUMMARY_FILE"

# Function to process each sample
process_sample() {
    local BASENAME=$1
    
    # Remove '_filtered' from BASENAME if it exists
    BASENAME=$(echo "$BASENAME" | sed 's/_filtered$//')
    # Check if files exist
    if [ ! -f "${RAW_READS_DIR}/${BASENAME}_1.fq.gz" ] || [ ! -f "${RAW_READS_DIR}/${BASENAME}_2.fq.gz" ]; then
        echo "Raw reads file ${RAW_READS_DIR}/${BASENAME}_1.fq.gz or ${RAW_READS_DIR}/${BASENAME}_2.fq.gz not found."
        return
    fi
    if [ ! -f "${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz" ] || [ ! -f "${FILTERED_DIR}/${BASENAME}_2_paired.fq.gz" ]; then
        echo "Filtered reads file ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz or ${FILTERED_DIR}/${BASENAME}_2_paired.fq.gz not found."
        return
    fi
    if [ ! -f "${FILTERED_DIR}/${BASENAME}_aligned_filtered.bam" ]; then
        echo "Filtered BAM file ${FILTERED_DIR}/${BASENAME}_aligned_filtered.bam not found."
        return
    fi
    if [ ! -f "${ANNOTATED_DIR}/${BASENAME}_aligned_filtered.emapper.annotations" ]; then
        echo "Annotation file ${ANNOTATED_DIR}/${BASENAME}_aligned_filtered.emapper.annotations not found."
        return
    fi

    # Calculate read counts
    RAW_READS=$(( ($(zcat "${RAW_READS_DIR}/${BASENAME}_1.fq.gz" | wc -l) + $(zcat "${RAW_READS_DIR}/${BASENAME}_2.fq.gz" | wc -l)) / 4 ))
    HQ_READS=$(( ($(zcat "${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz" | wc -l) + $(zcat "${FILTERED_DIR}/${BASENAME}_2_paired.fq.gz" | wc -l)) / 4 ))
    ALIGN_READS=$(samtools view -c "${FILTERED_DIR}/${BASENAME}_aligned_filtered.bam")
    ANNOT_READS=$(grep -c ">" "${ANNOTATED_DIR}/${BASENAME}_aligned_filtered.emapper.annotations")
    
    if [ "$RAW_READS" -eq 0 ]; then
        PROP_CONTAM="N/A"
        PROP_ANNOT="N/A"
    else
        PROP_CONTAM=$(echo "scale=2; 100 - ($ALIGN_READS / $RAW_READS * 100)" | bc)
        PROP_ANNOT=$(echo "scale=2; $ANNOT_READS / $RAW_READS * 100" | bc)
    fi
    
    # Append to summary file
    echo "$BASENAME,$RAW_READS,$HQ_READS,$ALIGN_READS,$ANNOT_READS,$PROP_CONTAM,$PROP_ANNOT" >> "$SUMMARY_FILE"
}

# Export the function and variables so they can be used by subshells
export -f process_sample
export RAW_READS_DIR FILTERED_DIR ANNOTATED_DIR SUMMARY_FILE

# Use parallel processing
ls "$RAW_READS_DIR"/*_1.fq.gz | while read file; do
    BASENAME=$(basename "$file" _1.fq.gz)
    ((i=i%NUM_THREADS)); ((i++==0)) && wait  # Limit the number of threads
    process_sample "$BASENAME" &  # Run each sample in the background
done

wait  # Wait for all background processes to finish before exiting the script

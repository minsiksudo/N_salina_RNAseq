#!/bin/bash

# Define the test sample (change this to the sample you want to test)
TEST_BASENAME="TN1506R0382"

# Define directories
RAW_READS_DIR="Nannochloropsis_20eax4_RNAseq"
FILTERED_DIR="filtered_reads"
ALIGN_DIR="aligned_reads"
COUNTS_DIR="counts"
QC_DIR="QC_results"
ADAPTERS_DIR="adapters"
REF_GENOME="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/GCA_004565275.1_ASM456527v1_genomic.fna"
HISAT2_INDEX="hisat2_index/nannochloropsis_salina"
GTF_FILE="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/genomic.gtf"
SUMMARY_FILE="results/read_counts_summary_test.csv"  # Define a test summary file for one sample
THREADS=40
JAVA_HEAP_SIZE="-Xms16g -Xmx64g"

# Create necessary output directories if they don't exist
mkdir -p $QC_DIR $FILTERED_DIR $ALIGN_DIR $COUNTS_DIR results

# Initialize summary file if it doesn't exist
if [ ! -f "$SUMMARY_FILE" ]; then
    echo "sample,raw_reads,high_quality_reads,forward_unpaired_reads,reverse_unpaired_reads,mapped_reads,assigned_reads" > "$SUMMARY_FILE"
fi

echo "Processing test sample $TEST_BASENAME..."

# Step 1: Run FastQC
echo "Running FastQC for $TEST_BASENAME..."
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate fastqc_env
fastqc -o $QC_DIR $RAW_READS_DIR/${TEST_BASENAME}_1.fq.gz $RAW_READS_DIR/${TEST_BASENAME}_2.fq.gz
conda deactivate

if [ $? -ne 0 ]; then
    echo "Error in FastQC step for $TEST_BASENAME"
    exit 1
fi

# Step 2: Run Trimmomatic
echo "Running Trimmomatic for $TEST_BASENAME..."
conda activate trimmomatic_env
java $JAVA_HEAP_SIZE -jar /home/bagel/miniconda3/envs/trimmomatic_env/share/trimmomatic-0.39-2/trimmomatic.jar PE -threads $THREADS \
    ${RAW_READS_DIR}/${TEST_BASENAME}_1.fq.gz ${RAW_READS_DIR}/${TEST_BASENAME}_2.fq.gz \
    ${FILTERED_DIR}/${TEST_BASENAME}_1_paired.fq.gz ${FILTERED_DIR}/${TEST_BASENAME}_1_unpaired.fq.gz \
    ${FILTERED_DIR}/${TEST_BASENAME}_2_paired.fq.gz ${FILTERED_DIR}/${TEST_BASENAME}_2_unpaired.fq.gz \
    ILLUMINACLIP:${ADAPTERS_DIR}/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
conda deactivate

if [ $? -ne 0 ]; then
    echo "Error in Trimmomatic step for $TEST_BASENAME"
    exit 1
fi

# Step 3: Run HISAT2
echo "Running HISAT2 for $TEST_BASENAME..."
conda activate hisat2_env
hisat2 -p $THREADS -x $HISAT2_INDEX \
    -1 ${FILTERED_DIR}/${TEST_BASENAME}_1_paired.fq.gz \
    -2 ${FILTERED_DIR}/${TEST_BASENAME}_2_paired.fq.gz \
    -S ${ALIGN_DIR}/${TEST_BASENAME}_aligned.sam
conda deactivate

if [ $? -ne 0 ]; then
    echo "Error in HISAT2 step for $TEST_BASENAME"
    exit 1
fi

# Step 4: Run Samtools Filtering, Sorting, and Indexing
echo "Running Samtools for $TEST_BASENAME..."
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate samtools_env

# Convert SAM to BAM, sort, and index
samtools view -@ $THREADS -bS ${ALIGN_DIR}/${TEST_BASENAME}_aligned.sam | \
samtools sort -@ $THREADS -o ${FILTERED_DIR}/${TEST_BASENAME}_sorted.bam
samtools index ${FILTERED_DIR}/${TEST_BASENAME}_sorted.bam

# Run Samtools flagstat to get alignment stats
echo "Generating alignment stats with Samtools flagstat for $TEST_BASENAME..."
SAMTOOLS_FLAGSTAT="${FILTERED_DIR}/${TEST_BASENAME}_flagstat.txt"
samtools flagstat ${FILTERED_DIR}/${TEST_BASENAME}_sorted.bam > $SAMTOOLS_FLAGSTAT

if [ ! -f "$SAMTOOLS_FLAGSTAT" ]; then
    echo "Error: Samtools flagstat file not found!"
    exit 1
fi

conda deactivate

# Step 5: Run featureCounts
echo "Running featureCounts for $TEST_BASENAME..."
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate subread_env
featureCounts -T $THREADS -p -a $GTF_FILE -o ${COUNTS_DIR}/${TEST_BASENAME}_counts.txt ${FILTERED_DIR}/${TEST_BASENAME}_sorted.bam

if [ $? -ne 0 ]; then
    echo "Error in featureCounts step for $TEST_BASENAME"
    exit 1
fi
conda deactivate

# Step 6: Generate summary statistics
echo "Generating summary for $TEST_BASENAME..."

# Count raw reads
RAW_READS=$(( ($(zcat ${RAW_READS_DIR}/${TEST_BASENAME}_1.fq.gz | wc -l) + $(zcat ${RAW_READS_DIR}/${TEST_BASENAME}_2.fq.gz | wc -l)) / 4 ))

# Count high-quality reads (after Trimmomatic)
HQ_READS=$(( ($(zcat ${FILTERED_DIR}/${TEST_BASENAME}_1_paired.fq.gz | wc -l) + $(zcat ${FILTERED_DIR}/${TEST_BASENAME}_2_paired.fq.gz | wc -l)) / 4 ))

# Extract forward-only and reverse-only unpaired reads from Trimmomatic output
FORWARD_ONLY_UNPAIRED=$(( $(zcat ${FILTERED_DIR}/${TEST_BASENAME}_1_unpaired.fq.gz | wc -l) / 4 ))
REVERSE_ONLY_UNPAIRED=$(( $(zcat ${FILTERED_DIR}/${TEST_BASENAME}_2_unpaired.fq.gz | wc -l) / 4 ))

# Extract mapped reads from flagstat
MAPPED_READS=$(grep "mapped (" ${SAMTOOLS_FLAGSTAT} | head -n 1 | cut -d " " -f 1)

# Extract assigned reads from featureCounts summary
FEATURECOUNTS_SUMMARY="${COUNTS_DIR}/${TEST_BASENAME}_counts.txt.summary"
ASSIGNED_READS=$(grep "Assigned" $FEATURECOUNTS_SUMMARY | awk '{print $2}')

# Append the stats to the summary file
echo "$TEST_BASENAME,$RAW_READS,$HQ_READS,$FORWARD_ONLY_UNPAIRED,$REVERSE_ONLY_UNPAIRED,$MAPPED_READS,$ASSIGNED_READS" >> "$SUMMARY_FILE"

echo "Pipeline completed successfully for $TEST_BASENAME!"



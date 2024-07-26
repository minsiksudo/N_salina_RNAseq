#!/bin/bash

# Initialize conda in the script
eval "$(conda shell.bash hook)"

# Directory paths
RAW_READS_DIR="Nannochloropsis_20eax4_RNAseq"
QC_DIR="results/qc_results"
FILTERED_DIR="filtered_reads"
HISAT2_INDEX_DIR="hisat2_index"
HISAT2_INDEX_PREFIX="reference_index"
ALIGN_DIR="aligned_reads"
ANNOTATED_DIR="annotated_reads"
REFERENCE_GENOME="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/GCA_004565275.1_ASM456527v1_genomic.fna"
REFERENCE_ANNOTATION="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/genomic.gff"
REFERENCE_GENES="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/cds_from_genomic.fna"
ADAPTER_FILE="adapters/TruSeq3-PE.fa" # Adjust this path
SUMMARY_FILE="results/read_counts_summary.csv"
NUM_THREADS=40  # Number of threads to use

# Create directories for results
mkdir -p ${QC_DIR}
mkdir -p ${FILTERED_DIR}
mkdir -p ${HISAT2_INDEX_DIR}
mkdir -p ${ALIGN_DIR}
mkdir -p ${ANNOTATED_DIR}

# Step 1: Quality Control (QC) with FastQC
echo "Running FastQC..."
conda activate base
fastqc -t ${NUM_THREADS} -o ${QC_DIR} ${RAW_READS_DIR}/*.fq.gz

# Step 2: Filter low-quality reads using Trimmomatic
echo "Filtering low-quality reads..."
conda activate trimmomatic_env
for READ1 in ${RAW_READS_DIR}/*_1.fq.gz; do
  READ2=${READ1/_1.fq.gz/_2.fq.gz}
  BASENAME=$(basename ${READ1} _1.fq.gz)
  
  # Trimming using Trimmomatic
  trimmomatic PE -threads ${NUM_THREADS} ${READ1} ${READ2} \
    ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz ${FILTERED_DIR}/${BASENAME}_1_unpaired.fq.gz \
    ${FILTERED_DIR}/${BASENAME}_2_paired.fq.gz ${FILTERED_DIR}/${BASENAME}_2_unpaired.fq.gz \
    ILLUMINACLIP:${ADAPTER_FILE}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
  
  # Counting the high-quality reads
  HIGH_QUALITY_READS=$(zcat ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz | echo $((`wc -l`/4)))
  
  # Update summary file
  echo "${BASENAME},${HIGH_QUALITY_READS}" >> ${SUMMARY_FILE}
done

# Step 3: Build the HISAT2 index
echo "Building HISAT2 index..."
conda activate base
hisat2-build ${REFERENCE_GENOME} ${HISAT2_INDEX_DIR}/${HISAT2_INDEX_PREFIX}

# Initialize summary file with headers
echo "Sample,Raw Reads,High Quality Reads,Reference Genome Aligned Reads,Annotated Reads,Proportion Contaminants,Proportion Annotated" > ${SUMMARY_FILE}

# Step 4: Align the reads with HISAT2 and process the alignments
for READ1 in ${FILTERED_DIR}/*_1_paired.fq.gz; do
  READ2=${READ1/_1_paired.fq.gz/_2_paired.fq.gz}
  BASENAME=$(basename ${READ1} _1_paired.fq.gz)
  
  echo "Processing ${BASENAME}..."
  
  # Align reads to the reference genome using HISAT2
  hisat2 -p ${NUM_THREADS} -x ${HISAT2_INDEX_DIR}/${HISAT2_INDEX_PREFIX} -1 ${READ1} -2 ${READ2} -S ${ALIGN_DIR}/${BASENAME}_aligned.sam
  
  # Convert SAM to BAM, sort, and index
  samtools view -S -b ${ALIGN_DIR}/${BASENAME}_aligned.sam > ${ALIGN_DIR}/${BASENAME}_aligned.bam
  samtools sort ${ALIGN_DIR}/${BASENAME}_aligned.bam -o ${ALIGN_DIR}/${BASENAME}_aligned_sorted.bam
  samtools index ${ALIGN_DIR}/${BASENAME}_aligned_sorted.bam
  
  # Extract aligned reads to the filtered directory
  samtools view -b -F 4 ${ALIGN_DIR}/${BASENAME}_aligned_sorted.bam > ${FILTERED_DIR}/${BASENAME}_aligned_filtered.bam
  
  # Count raw reads, high-quality reads, and reference genome aligned reads
  RAW_READS=$(zcat ${RAW_READS_DIR}/${BASENAME}_1.fq.gz | echo $((`wc -l`/4)))
  HIGH_QUALITY_READS=$(zcat ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz | echo $((`wc -l`/4)))
  REFERENCE_ALIGNED_READS=$(samtools view -c -F 4 ${ALIGN_DIR}/${BASENAME}_aligned_sorted.bam)
  
  # Update summary file with raw reads, high-quality reads, and reference genome aligned reads
  sed -i "s/^${BASENAME},.*/&${RAW_READS},${HIGH_QUALITY_READS},${REFERENCE_ALIGNED_READS}/" ${SUMMARY_FILE}
done

# Step 5: Quantify transcripts using featureCounts
echo "Running featureCounts..."
featureCounts -T ${NUM_THREADS} -a ${REFERENCE_ANNOTATION} -o ${ANNOTATED_DIR}/counts.txt ${FILTERED_DIR}/*_aligned_filtered.bam

# Step 6: Functional annotation with eggNOG-mapper
echo "Running eggNOG-mapper..."
conda activate eggnog_mapper_env
emapper.py -i ${REFERENCE_GENES} -o ${ANNOTATED_DIR}/annotations --cpu ${NUM_THREADS}

# Calculate proportions and update summary file
for READ1 in ${FILTERED_DIR}/*_1_paired.fq.gz; do
  BASENAME=$(basename ${READ1} _1_paired.fq.gz)
  
  ANNOTATED_READS=$(grep -c "^" ${ANNOTATED_DIR}/counts.txt)
  PROPORTION_CONTAMINANTS=$(awk "BEGIN {printf \"%.2f\", (100 - (${REFERENCE_ALIGNED_READS}/${RAW_READS})*100)}")
  PROPORTION_ANNOTATED=$(awk "BEGIN {printf \"%.2f\", (${ANNOTATED_READS}/${REFERENCE_ALIGNED_READS})*100}")
  
  # Update summary file
  sed -i "s/^${BASENAME},.*/&${ANNOTATED_READS},${PROPORTION_CONTAMINANTS},${PROPORTION_ANNOTATED}/" ${SUMMARY_FILE}
done

echo "Alignment, filtering, and annotation completed. Summary file generated: ${SUMMARY_FILE}"

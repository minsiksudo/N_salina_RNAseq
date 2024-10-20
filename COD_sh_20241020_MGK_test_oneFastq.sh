#!/bin/bash

# Define the sample to test with (replace 'test_sample' with your sample name)
TEST_SAMPLE="TN1511R0955--CCGTCCCG"

# Define directories
RAW_READS_DIR="Nannochloropsis_20eax4_RNAseq"
FILTERED_DIR="filtered_reads"
ALIGN_DIR="aligned_reads"
ANNOTATED_DIR="annotated_reads"
QC_DIR="QC_results"
COUNTS_DIR="counts"
FASTA_DIR="fasta_reads"

# Create necessary directories if they don't exist
mkdir -p $FILTERED_DIR $ALIGN_DIR $ANNOTATED_DIR $QC_DIR $COUNTS_DIR $FASTA_DIR

# Check if the necessary environments are installed and activated
echo "Checking conda environments..."

# FastQC step
echo "Running FastQC on $TEST_SAMPLE..."
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate fastqc_env
fastqc -o $QC_DIR ${RAW_READS_DIR}/${TEST_SAMPLE}_1.fq.gz ${RAW_READS_DIR}/${TEST_SAMPLE}_2.fq.gz
if [ $? -ne 0 ]; then
    echo "FastQC step failed"
    exit 1
fi
conda deactivate

# Trimmomatic step
echo "Running Trimmomatic on $TEST_SAMPLE..."
conda activate trimmomatic_env
java -Xms16g -Xmx64g -jar /home/bagel/miniconda3/envs/trimmomatic_env/share/trimmomatic-0.39-2/trimmomatic.jar PE \
    ${RAW_READS_DIR}/${TEST_SAMPLE}_1.fq.gz ${RAW_READS_DIR}/${TEST_SAMPLE}_2.fq.gz \
    ${FILTERED_DIR}/${TEST_SAMPLE}_1_paired.fq.gz ${FILTERED_DIR}/${TEST_SAMPLE}_1_unpaired.fq.gz \
    ${FILTERED_DIR}/${TEST_SAMPLE}_2_paired.fq.gz ${FILTERED_DIR}/${TEST_SAMPLE}_2_unpaired.fq.gz \
    ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
if [ $? -ne 0 ]; then
    echo "Trimmomatic step failed"
    exit 1
fi
conda deactivate

# HISAT2 step
echo "Running HISAT2 on $TEST_SAMPLE..."
conda activate hisat2_env
hisat2 -p 20 -x hisat2_index/nannochloropsis_salina -1 ${FILTERED_DIR}/${TEST_SAMPLE}_1_paired.fq.gz -2 ${FILTERED_DIR}/${TEST_SAMPLE}_2_paired.fq.gz -S ${ALIGN_DIR}/${TEST_SAMPLE}_aligned.sam
if [ $? -ne 0 ]; then
    echo "HISAT2 step failed"
    exit 1
fi
conda deactivate

# Samtools filtering step
echo "Running Samtools filtering on $TEST_SAMPLE..."
conda activate samtools_env
samtools view -@ 40 -bS ${ALIGN_DIR}/${TEST_SAMPLE}_aligned.sam | samtools sort -@ 40 -m 8G -o ${FILTERED_DIR}/${TEST_SAMPLE}_sorted.bam
if [ $? -ne 0 ]; then
    echo "Samtools sort failed"
    exit 1
fi
samtools index ${FILTERED_DIR}/${TEST_SAMPLE}_sorted.bam
if [ $? -ne 0 ]; then
    echo "Samtools index failed"
    exit 1
fi
samtools view -@ 40 -b -q 20 ${FILTERED_DIR}/${TEST_SAMPLE}_sorted.bam > ${FILTERED_DIR}/${TEST_SAMPLE}_filtered.bam
if [ $? -ne 0 ]; then
    echo "Samtools filtering failed"
    exit 1
fi
samtools index ${FILTERED_DIR}/${TEST_SAMPLE}_filtered.bam
if [ $? -ne 0 ]; then
    echo "Samtools filtering failed"
    exit 1
fi
conda deactivate

# EggNOG-mapper step
echo "Running EggNOG-mapper on $TEST_SAMPLE..."
conda activate eggnog_mapper_env
samtools fasta ${FILTERED_DIR}/${TEST_SAMPLE}_filtered.bam > ${FASTA_DIR}/${TEST_SAMPLE}.fasta
if [ $? -ne 0 ]; then
    echo "Conversion to FASTA failed"
    exit 1
fi
emapper.py --cpu 40 --output ${ANNOTATED_DIR}/${TEST_SAMPLE} --data_dir /mnt/4T_samsung/Dropbox/Sequencing_archive/eggnog_db -i ${FASTA_DIR}/${TEST_SAMPLE}.fasta --override --translate -m diamond
if [ $? -ne 0 ]; then
    echo "EggNOG-mapper step failed"
    exit 1
fi
conda deactivate

echo "Test pipeline for $TEST_SAMPLE completed successfully!"

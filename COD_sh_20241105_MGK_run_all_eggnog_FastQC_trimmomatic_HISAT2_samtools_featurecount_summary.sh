#!/bin/bash

# Define directories
RAW_READS_DIR="Nannochloropsis_20eax4_RNAseq"
FILTERED_DIR="filtered_reads"
ALIGN_DIR="aligned_reads"
COUNTS_DIR="eggnog_counts"
QC_DIR="QC_results"
ADAPTERS_DIR="adapters"
REF_GENOME="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/GCA_004565275.1_ASM456527v1_genomic.fna"
HISAT2_INDEX="hisat2_index/nannochloropsis_salina"  # Corrected basename
GFF_FILE="eggnog_output/nannochloropsis_with_eggnog.gtf"
SUMMARY_FILE="results/eggnog_read_counts_summary.csv"  # Define the summary file
THREADS=40
JAVA_HEAP_SIZE="-Xms16g -Xmx64g"

# Create necessary output directories if they don't exist
mkdir -p $QC_DIR $FILTERED_DIR $ALIGN_DIR $COUNTS_DIR results

# Initialize summary file if it doesn't exist
if [ ! -f "$SUMMARY_FILE" ]; then
    echo "sample,raw_reads,high_quality_reads,forward_unpaired_reads,reverse_unpaired_reads,mapped_reads,assigned_reads" > "$SUMMARY_FILE"
fi

# Fix GTF file
echo "Fixing GTF file..."
awk 'BEGIN{FS=OFS="\t"} 
{
    gsub(/;;/, ";", $9)
    gsub(/;  gene_id/, "; gene_id", $9)
    if ($9 !~ /gene_id/) {
        split($9, arr, "transcript_id ")
        split(arr[2], arr2, ";")
        transcript_id = arr2[1]
        $9 = $9 "; gene_id " transcript_id
    }
    print
}' $GFF_FILE > ${GFF_FILE%.gtf}_fixed.gtf

# Clean up GTF file
sed 's/;;/;/g' ${GFF_FILE%.gtf}_fixed.gtf > ${GFF_FILE%.gtf}_fixed_clean.gtf

FEATURECOUNTS_INPUT=${GFF_FILE%.gtf}_fixed_clean.gtf

# Loop through all paired FASTQ files in the raw reads directory
for FILE in $RAW_READS_DIR/*_1.fq.gz; do
    # Extract the sample name (basename without _1.fq.gz or _2.fq.gz)
    BASENAME=$(basename $FILE _1.fq.gz)
    
    echo "Processing sample $BASENAME..."

    # Step 1: Run FastQC
    echo "Running FastQC for $BASENAME..."
    source /home/bagel/miniconda3/etc/profile.d/conda.sh
    conda activate fastqc_env
    fastqc -o $QC_DIR $RAW_READS_DIR/${BASENAME}_1.fq.gz $RAW_READS_DIR/${BASENAME}_2.fq.gz
    conda deactivate

    if [ $? -ne 0 ]; then
        echo "Error in FastQC step for $BASENAME"
        exit 1
    fi

    # Step 2: Run Trimmomatic
    echo "Running Trimmomatic for $BASENAME..."
    conda activate trimmomatic_env
    java $JAVA_HEAP_SIZE -jar /home/bagel/miniconda3/envs/trimmomatic_env/share/trimmomatic-0.39-2/trimmomatic.jar PE -threads $THREADS \
        ${RAW_READS_DIR}/${BASENAME}_1.fq.gz ${RAW_READS_DIR}/${BASENAME}_2.fq.gz \
        ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz ${FILTERED_DIR}/${BASENAME}_1_unpaired.fq.gz \
        ${FILTERED_DIR}/${BASENAME}_2_paired.fq.gz ${FILTERED_DIR}/${BASENAME}_2_unpaired.fq.gz \
        ILLUMINACLIP:${ADAPTERS_DIR}/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
    conda deactivate

    if [ $? -ne 0 ]; then
        echo "Error in Trimmomatic step for $BASENAME"
        exit 1
    fi

    # Step 3: Run HISAT2 (corrected index basename)
    echo "Running HISAT2 for $BASENAME..."
    conda activate hisat2_env
    hisat2 -p $THREADS -x $HISAT2_INDEX \
        -1 ${FILTERED_DIR}/${BASENAME}_1_paired.fq.gz \
        -2 ${FILTERED_DIR}/${BASENAME}_2_paired.fq.gz \
        -S ${ALIGN_DIR}/${BASENAME}_aligned.sam
    conda deactivate

    if [ $? -ne 0 ]; then
        echo "Error in HISAT2 step for $BASENAME"
        exit 1
    fi

    # Step 4: Run Samtools Filtering, Sorting, and Indexing (correct environment setup)
    echo "Running Samtools for $BASENAME..."
    source /home/bagel/miniconda3/etc/profile.d/conda.sh
    conda activate samtools_env

    samtools view -@ $THREADS -bS ${ALIGN_DIR}/${BASENAME}_aligned.sam | \
    samtools sort -@ $THREADS -o ${FILTERED_DIR}/${BASENAME}_sorted.bam
    samtools index ${FILTERED_DIR}/${BASENAME}_sorted.bam

    SAMTOOLS_FLAGSTAT="${FILTERED_DIR}/${BASENAME}_flagstat.txt"
    samtools flagstat ${FILTERED_DIR}/${BASENAME}_sorted.bam > $SAMTOOLS_FLAGSTAT

    if [ ! -f "$SAMTOOLS_FLAGSTAT" ]; then
        echo "Error: Samtools flagstat file not found!"
        exit 1
    fi

    conda deactivate

    # Step 5: Run featureCounts (adjusted for paired-end data)
    echo "Running featureCounts for $BASENAME..."
    conda activate subread_env
    
    featureCounts -T $THREADS -p -B -C \
                  -a "$FEATURECOUNTS_INPUT" \
                  -o ${COUNTS_DIR}/${BASENAME}_counts.txt \
                  ${FILTERED_DIR}/${BASENAME}_sorted.bam

    if [ $? -ne 0 ]; then
        echo "Error in featureCounts step for $BASENAME"
        exit 1
    fi
    
    conda deactivate
    
done

echo "Pipeline completed successfully!"


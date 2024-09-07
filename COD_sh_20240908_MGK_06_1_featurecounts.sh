# Example script to run featureCounts on BAM files
GTF_FILE="path/to/Nannochloropsis_salina.gtf"
OUTPUT_DIR="counts"
THREADS=4
BAM_DIR="aligned_reads"

mkdir -p $OUTPUT_DIR

for BAM in $BAM_DIR/*.bam; do
    BASENAME=$(basename $BAM .bam)
    featureCounts -a $GTF_FILE -o ${OUTPUT_DIR}/${BASENAME}_counts.txt -T $THREADS $BAM
done


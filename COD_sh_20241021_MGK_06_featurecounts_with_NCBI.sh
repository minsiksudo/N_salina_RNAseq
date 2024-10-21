# COD_sh_20240908_MGK_06_featurecounts_with_NCBI.sh

# Define directories
ALIGN_DIR="aligned_reads"
COUNTS_DIR="counts"
GTF_FILE="ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/annotations.gtf"  # Update with your GTF file path
THREADS=20  # Adjust based on your system

# Create output directory if it doesn't exist
mkdir -p $COUNTS_DIR

# Activate the environment with featureCounts (if required)
source /home/bagel/miniconda3/etc/profile.d/conda.sh
conda activate subread_env

# Run featureCounts on all BAM files
for bamfile in $ALIGN_DIR/*_sorted.bam; do
    BASENAME=$(basename $bamfile _sorted.bam)
    echo "Processing $BASENAME..."
    featureCounts -T $THREADS -a $GTF_FILE -o ${COUNTS_DIR}/${BASENAME}_counts.txt $bamfile
    if [ $? -ne 0 ]; then
        echo "Error in featureCounts for $BASENAME"
        continue
    fi
done

# Deactivate the environment
conda deactivate

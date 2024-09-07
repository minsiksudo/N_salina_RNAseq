# Set directories
COUNTS_DIR="counts"
ANNOTATED_DIR="annotated_reads"
OUTPUT_DIR="combined_outputs"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through each counts file
for file in $COUNTS_DIR/*_counts.txt; do
    BASENAME=$(basename $file _counts.txt)
    
    echo "Processing $BASENAME..."
    
    # Path to EggNOG annotation file
    ANNOTATION_FILE="$ANNOTATED_DIR/${BASENAME}_filtered.emapper.annotations"
    
    if [ -f "$ANNOTATION_FILE" ]; then
        echo "Found annotation for $BASENAME. Combining read counts with EggNOG annotations..."
        
        # Join counts and EggNOG annotations by gene (assuming gene_id is in the first column of both)
        awk 'BEGIN {FS=OFS="\t"} FNR==NR {a[$1]=$0; next} $1 in a {print a[$1], $0}' $ANNOTATION_FILE $file > ${OUTPUT_DIR}/${BASENAME}_combined.txt
    else
        echo "Warning: Annotation file for $BASENAME not found. Skipping..."
    fi
done

echo "All files processed. Combined outputs are in the $OUTPUT_DIR directory."

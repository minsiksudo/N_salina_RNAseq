#First modification  to make gene_id compativle with featurecounts

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
}' nannochloropsis_annotation.gtf > nannochloropsis_annotation_fixed2.gtf

#Second modification

sed 's/;;/;/g' nannochloropsis_annotation_fixed2.gtf > nannochloropsis_annotation_fixed3.gtf

# Testing one file to see if it works
featureCounts -T 40 -p -B -C -a nannochloropsis_annotation_fixed3.gtf \
              -o TN1506R0382_counts.txt aligned_reads/TN1506R0382_aligned_sorted.bam


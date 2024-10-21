import csv

# Paths to input and output files
eggnog_annotations_file = "eggnog_output/nannochloropsis_annotation.emapper.annotations"
original_gtf_file = "ref_genome/GCA_004565275.1/ncbi_dataset/data/GCA_004565275.1/genomic.gtf"
output_gtf_file = "nannochloropsis_with_eggnog.gtf"

# Dictionary to store functional annotations (EggNOG)
annotations = {}

# Read EggNOG annotations and store in dictionary
with open(eggnog_annotations_file, 'r') as eggnog_file:
    reader = csv.reader(eggnog_file, delimiter='\t')
    for row in reader:
        # Skip header lines
        if row[0].startswith("#") or len(row) < 6:
            continue
        query_id = row[0]
        description = row[5]  # Functional description
        annotations[query_id] = description

# Open the original GTF and output GTF file
with open(original_gtf_file, 'r') as gtf_in, open(output_gtf_file, 'w') as gtf_out:
    for line in gtf_in:
        # If the line is a comment or metadata, just copy it
        if line.startswith("#"):
            gtf_out.write(line)
            continue
        
        # Parse the GTF fields
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue  # Skip malformed lines
        
        # Extract gene or transcript ID from attributes field
        attributes = fields[8]
        gene_id = None
        if "gene_id" in attributes:
            gene_id = attributes.split("gene_id")[1].split(";")[0].replace('"', '').strip()

        # Add EggNOG functional annotations if available
        if gene_id and gene_id in annotations:
            description = annotations[gene_id]
            new_attributes = f'{attributes} functional_description "{description}";'
            fields[8] = new_attributes
        
        # Write the updated line to the output GTF
        gtf_out.write("\t".join(fields) + "\n")

print(f"New GTF file with EggNOG annotations saved to {output_gtf_file}")

**QCing, mapping, and annotating RNAs of** ***Nannochloropsis salina*** **CCMP1776** were conducted for a collaborative work with **Prof. Nam Kyu Kang at KHU**. Plus, Minsik's lab is encharge of statistical analyses to find out how the microalgal pathways expressed differently based on different culture conditions, by time. 

**Materials and Methods**

*Quality Control and Preprocessing*

Quality control, trimming and filtering, alignment, annotation and quantification

The quality of raw RNA-seq reads was assessed using FastQC v0.11.9 (Andrews, 2010). Reads were evaluated for various quality metrics including per base sequence quality, GC content, and adapter contamination. Low-quality reads and adapter sequences were removed using Trimmomatic v0.39 (Bolger et al., 2014). The following parameters were used: ILLUMINACLIP:TruSeq3-PE.fa:2:30:10, SLIDINGWINDOW:4:20, and MINLEN:36. High-quality paired-end reads were retained for downstream analysis. The Nannochloropsis salina CCMP1776 reference genome (GCA_004565275.1) and corresponding annotation files were obtained from the NCBI database. The genome was indexed using HISAT2 v2.2.1 (Kim et al., 2015). High-quality reads were aligned to the reference genome using HISAT2 with the following parameters: -p 48 for 48 threads, ensuring efficient utilization of computational resources. Aligned reads were converted to BAM format, sorted, and indexed using Samtools v1.10 (Li et al., 2009). Gene-level read counts were generated using featureCounts v2.0.1 (Liao et al., 2014) with the reference annotation file. The parameters -T 48 were used to specify the number of threads. Functional annotation of the gene models was performed using eggNOG-mapper v2 (Huerta-Cepas et al., 2017). Annotated reads were then quantified, and the proportion of contaminants was determined based on the alignment statistics.

Statistical Analysis

…TBA

Data Availability

The raw sequencing data and processed files are available in the NCBI Sequence Read Archive under the accession number [SRA OOOOOOO (TBA); PRJN OOOOOOO (TBA)]. Scripts used for the analysis are available at [https://github.com/minsiksudo/N_salina_RNAseq].

References:

Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. Available online: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.

Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. Nature Methods, 12(4), 357-360.

Li, H., Handsaker, B., Wysoker, A., et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079.

Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930.

Huerta-Cepas, J., Forslund, K., Szklarczyk, D., et al. (2017). eggNOG 4.5: a hierarchical orthology framework with improved functional annotations for eukaryotic, prokaryotic and viral sequences. Nucleic Acids Research, 45(D1), D309-D314.



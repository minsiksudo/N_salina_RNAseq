**QCing, mapping, and annotating RNAs of** ***Nannochloropsis salina*** **CCMP1776** were conducted for a collaborative work with **Prof. Nam Kyu Kang at KHU**. Plus, Minsik's lab is encharge of statistical analyses to find out how the microalgal pathways expressed differently based on different culture conditions, by time. 

**Materials and Methods**

*Quality control and data reprocessing*

Quality control, trimming, alignment, filtering, annotation and quantification

The quality of raw RNA-seq reads was initially assessed using FastQC v0.11.9 [^1] to evaluate metrics such as per base sequence quality, GC content, and adapter contamination. Low-quality reads and adapter sequences were removed using Trimmomatic v0.39 [^2] with parameters: ILLUMINACLIP:TruSeq3-PE.fa:2:30:10, SLIDINGWINDOW:4:20, and MINLEN:36. High-quality paired-end reads were retained for downstream analysis. The Nannochloropsis salina CCMP1776 reference genome (GCA_004565275.1) and the corresponding annotation files were obtained from NCBI. The genome was indexed using HISAT2 v2.2.1 [^3] with [include any specific parameters]. High-quality reads were aligned to the reference genome using HISAT2 with the parameter -p 48 for efficient multi-threading. Aligned reads were converted to BAM format, sorted, and indexed using Samtools v1.10 [^4]. Gene-level quantification was performed using featureCounts v2.0.1 [^5] with parameters -T 48 to specify the number of threads. Functional annotation of gene models was performed using eggNOG-mapper v2 [^6] with database version 5.0.0. Annotated reads were then quantified, and the proportion of contaminants was determined based on alignment statistics and post-alignment QC results. Post-alignment QC was performed using tools such as FastQC, Picard [^7], SAMtools, and Qualimap [^8] to ensure the quality of alignments and identify potential issues.


Statistical Analysis

…TBA

Data Availability

The raw sequencing data and processed files are available in the NCBI Sequence Read Archive under the accession number [SRA OOOOOOO (TBA); PRJN OOOOOOO (TBA)]. Scripts used for the analysis are available at [https://github.com/minsiksudo/N_salina_RNAseq].

References

	1.	Andrews, S. (2010). FastQC: A Quality Control tool for High Throughput Sequence Data. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.
	2.	Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120. doi:10.1093/bioinformatics/btu170.
	3.	Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: A fast spliced aligner with low memory requirements. Nature Methods, 12(4), 357-360. doi:10.1038/nmeth.3317.
	4.	Li, H., Handsaker, B., Wysoker, A., et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079. doi:10.1093/bioinformatics/btp352.
	5.	Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: An efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930. doi:10.1093/bioinformatics/btt656.




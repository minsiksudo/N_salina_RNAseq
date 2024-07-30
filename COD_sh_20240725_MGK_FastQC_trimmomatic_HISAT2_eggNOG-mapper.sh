#!/bin/bash

# Run FastQC
#echo "Running FastQC..."
#./COD_sh_20240726_MGK_01_fastqc.sh
#if [ $? -ne 0 ]; then
#    echo "Error in FastQC step"
#    exit 1
#fi

# Run Trimmomatic
#echo "Running Trimmomatic..."
#./COD_sh_20240726_MGK_02_trimmomatic.sh
#if [ $? -ne 0 ]; then
#    echo "Error in Trimmomatic step"
#    exit 1
#fi

# Run HISAT2
#echo "Running HISAT2..."
#./COD_sh_20240726_MGK_03_HISAT2.sh
#if [ $? -ne 0 ]; then
#    echo "Error in HISAT2 step"
#    exit 1
#fi

# Run Samtools Filtering
#echo "Running Samtools Filtering..."
#./COD_sh_20240726_MGK_04_sam_tools_filtering.sh
#if [ $? -ne 0 ]; then
#    echo "Error in Samtools Filtering step"
#    exit 1
#fi

# Run EggNOG-mapper
echo "Running EggNOG-mapper..."
./COD_sh_20240726_MGK_05_eggnog_mapper.sh
if [ $? -ne 0 ]; then
    echo "Error in EggNOG-mapper step"
    exit 1
fi

# Run featureCounts

echo "Running fetureCounts..."
./COD_sh_20240730_MGK_06_featurecounts.sh
if [ $? -ne 0 ]; then
    echo "Error in featureCount step"
    exit 1
fi



# Generate Summary
echo "Generating Summary..."
./COD_sh_20240726_MGK_06_summary.sh
if [ $? -ne 0 ]; then
    echo "Error in Summary generation step"
    exit 1
fi

echo "All steps completed successfully!"

#!/bin/bash

# Step 1: Run FastQC
echo "Running FastQC..."
./COD_sh_20240726_MGK_01_fastqc.sh
if [ $? -ne 0 ]; then
    echo "Error in FastQC step"
    exit 1
fi

# Step 2: Run Trimmomatic
echo "Running Trimmomatic..."
./COD_sh_20240726_MGK_02_trimmomatic.sh
if [ $? -ne 0 ]; then
    echo "Error in Trimmomatic step"
    exit 1
fi

# Step 3: Run HISAT2
echo "Running HISAT2..."
./COD_sh_20240726_MGK_03_HISAT2.sh
if [ $? -ne 0 ]; then
    echo "Error in HISAT2 step"
    exit 1
fi

# Step 4: Run Samtools Filtering
echo "Running Samtools Filtering..."
./COD_sh_20240726_MGK_04_sam_tools_filtering.sh
if [ $? -ne 0 ]; then
    echo "Error in Samtools Filtering step"
    exit 1
fi

# Step 5: Run featureCounts using NCBI Annotation
echo "Running featureCounts..."
./COD_sh_20240908_MGK_06_featurecounts_with_NCBI.sh
if [ $? -ne 0 ]; then
    echo "Error in featureCounts step"
    exit 1
fi

# Step 6: Generate Summary
echo "Generating Summary..."
./COD_sh_20240726_MGK_07_summary.sh
if [ $? -ne 0 ]; then
    echo "Error in Summary generation step"
    exit 1
fi

echo "All steps completed successfully!"

#!/bin/bash
# IEKNGHER1ujm

# Setup
module load bioinfo/bcftools
DUMP=/home/dwragg/work/Analysis

# VCF arrays
readarray VCFS_A < ${DUMP}/vcfs/lists/A_vcfs.list
readarray VCFS_M < ${DUMP}/vcfs/lists/M_vcfs.list
readarray VCFS_C < ${DUMP}/vcfs/lists/C_vcfs.list
readarray VCFS_O < ${DUMP}/vcfs/lists/O_vcfs.list
readarray VCFS_JFM < ${DUMP}/vcfs/lists/JFM_vcfs.list
readarray VCFS_OTH < ${DUMP}/vcfs/lists/OTH_vcfs.list
tmp=(${VCFS_A[@]} ${VCFS_M[@]} ${VCFS_C[@]} ${VCFS_O[@]} ${VCFS_JFM[@]} ${VCFS_OTH[@]})

# Convert array to string
AMCO=${tmp[@]}


# Identify intersecting sites across all core VCF files
bcftools isec -n+${#tmp[@]} -O z -w1 -o ${DUMP}/vcfs/test.vcf.gz ${AMCO}




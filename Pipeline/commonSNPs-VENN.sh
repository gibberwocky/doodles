#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

# Paths to main folders
DUMP="/home/dwragg/work/Analysis"
PIPE="/home/dwragg/work/Pipeline"
PERL="/save/dwragg/Scripts/Perl"
RSCRIPTS="/save/seqapipop/Scripts/R"

# The aim of this script is to identify SNPs in each lienage/population that
# are unique/shared by means of Venn diagram. Note, this needs doing prior
# to imputation otherwise everyone will have the same number of SNPs!

# Specify set of master sites
#sites=${DUMP}/sites/HARP-d9-3113.vcf.gz
#POP="A"
sites=${DUMP}/sites/JFM-OTH-d9-1278.vcf.gz
POP="JFMOTH"

bcftools merge -l ${DUMP}/vcfs/lists/${POP}_vcfs.list \
  | bcftools view --types snps \
  -T ${sites} \
  -M 2 -O v -o ${DUMP}/vcfs/${POP}_raw.vcf

#			n
# A_raw.vcf		8007044
# M_raw.vcf		2635813
# C_raw.vcf		2277921	
# O_raw.vcf		3204301
# JFM_raw.vcf		3417776
# OTH_raw.vcf		2846383	
# JFMOTH_raw.vcf	3684050

# Generate stats and Venn Diagram (max 5 VCF files possible)
SAMPLE="raw"
${RSCRIPTS}/VennDiagram-raw.R ${SAMPLE} ${DUMP}/vcfs/HARP-JFM-OTH-07112014/${SAMPLE}_VCF_VENN.tiff \
  ${DUMP}/vcfs/HARP-JFM-OTH-07112014/A_raw.vcf \
  ${DUMP}/vcfs/HARP-JFM-OTH-07112014/M_raw.vcf \
  ${DUMP}/vcfs/HARP-JFM-OTH-07112014/C_raw.vcf \
  ${DUMP}/vcfs/HARP-JFM-OTH-07112014/O_raw.vcf \
  ${DUMP}/vcfs/HARP-JFM-OTH-07112014/JFMOTH_raw.vcf 






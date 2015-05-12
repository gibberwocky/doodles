#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

# Paths to main folders
DUMP="/home/dwragg/work/Analysis"
PIPE="/home/dwragg/work/Pipeline"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"

# AMEL-4M
OUT=${DUMP}/vcfs/HARP-JFM-OTH-07112014


# Merge data (AMCO and haploids, for checking imputation rate)
# Data pre-imputation (flash format first due to PL field error)
bcftools annotate -x format -O z -o ${OUT}/A_xformat.vcf.gz ${OUT}/A_subhet.vcf.gz
bcftools annotate -x format -O z -o ${OUT}/M_xformat.vcf.gz ${OUT}/M_subhet.vcf.gz
bcftools annotate -x format -O z -o ${OUT}/C_xformat.vcf.gz ${OUT}/C_subhet.vcf.gz
bcftools annotate -x format -O z -o ${OUT}/O_xformat.vcf.gz ${OUT}/O_subhet.vcf.gz
bcftools annotate -x format -O z -o ${OUT}/JFM_xformat.vcf.gz ${OUT}/JFM_subhet.vcf.gz
bcftools annotate -x format -O z -o ${OUT}/OTH_xformat.vcf.gz ${OUT}/OTH_subhet.vcf.gz
tabix -fp vcf ${OUT}/A_xformat.vcf.gz
tabix -fp vcf ${OUT}/M_xformat.vcf.gz
tabix -fp vcf ${OUT}/C_xformat.vcf.gz
tabix -fp vcf ${OUT}/O_xformat.vcf.gz
tabix -fp vcf ${OUT}/JFM_xformat.vcf.gz
tabix -fp vcf ${OUT}/OTH_xformat.vcf.gz
bcftools merge -m snps -O z \
  ${OUT}/A_xformat.vcf.gz \
  ${OUT}/M_xformat.vcf.gz \
  ${OUT}/C_xformat.vcf.gz \
  ${OUT}/O_xformat.vcf.gz \
  ${OUT}/JFM_xformat.vcf.gz \
  ${OUT}/OTH_xformat.vcf.gz \
  -o ${OUT}/S2_AMCO-hapJFMOTH_xformat.vcf.gz
tabix -fp vcf ${OUT}/S2_AMCO-hapJFMOTH_xformat.vcf.gz

# Data post-imputation
bcftools merge -m snps -O z \
  ${OUT}/A_beagle.vcf.gz \
  ${OUT}/M_beagle.vcf.gz \
  ${OUT}/C_beagle.vcf.gz \
  ${OUT}/O_beagle.vcf.gz \
  ${OUT}/JFM_beagle.vcf.gz \
  ${OUT}/OTH_beagle.vcf.gz \
  -o ${OUT}/S2_AMCO-hapJFMOTH_beagle.vcf.gz
tabix -fp vcf ${OUT}/S2_AMCO-hapJFMOTH_beagle.vcf.gz

# Replace underscores in sample IDs
FILE=S2_AMCO-hapJFMOTH_beagle.vcf
gunzip ${OUT}/${FILE}.gz
sed -i '/^#CHROM/s/_/-/g' ${OUT}/${FILE}

# Check genotyping rate
${PLINK}/plink --vcf ${OUT}/${FILE} \
  --allow-no-sex \
  --allow-extra-chr \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${OUT}/${FILE} \
  --missing




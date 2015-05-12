#!/bin/bash
# IEKNGHER1ujm

# Set-up
module load bioinfo/bcftools
module load bioinfo/Java7 
BEDT="/usr/local/bioinfo/src/bedtools/bedtools2-2.20.1/bin"
VEP="/usr/local/bioinfo/src/ensembl-api/current/ensembl-tools-release-78/scripts/variant_effect_predictor"

# Data
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
SITES=${DUMP}/vcfs/HARP-JFM-OTH-15012015
VCF=JFMOTH_preImp_chrs16_qc.vcf.gz


# Results
BEDS=(${DUMP}/selection/results/raw/*.bed)

# Clean BED file by filtering out solo hits (likely false-positives)
for ibed in ${BEDS[@]}
do
  ${BEDT}/mergeBed -d 1000 -c 1,5,5,5 -o count,distinct,min,max \
    -i ${ibed} | awk '{ if ($4 > 1) print $0 }' - > ${ibed}.clean
  bcftools view --types snps -M 2 -O v \
    -R ${ibed}.clean \
    -o ${ibed}.clean.vcf \
    ${SITES}/${VCF}
done

# Except for fst.bed which can support solo hit windows (fewer hits)
ibed="fst.bed"
${BEDT}/mergeBed -d 1000 -c 1,5,5,5 -o count,distinct,min,max \
  -i ${ibed} > ${ibed}.clean
bcftools view --types snps -M 2 -O v \
  -R ${ibed}.clean \
  -o ${ibed}.clean.vcf \
  ${SITES}/${VCF}


# Run VEP on the resulting VCF files
VCFS=(${DUMP}/selection/results/raw/*.vcf)
mkdir -p ${DUMP}/logs/VEP-selection
for ivcf in ${VCFS[@]}
do
  qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/VEP-selection \
    -e ${DUMP}/logs/VEP-selection \
    ${PIPE}/vep.sh \
      -i ${ivcf} \
      -o ${ivcf}
done






# 1) Import VEP outputs into Calc
# 2) Copy Genes, paste to new sheet and filter (no field, remove duplicates)
# 3) Use gene list as input to Biomart to recover:
#	GO terms (accession, name, domain)
#	Homology (drosophila gene ID, homology type, &identity, &dmel identity, confidence)
# 4) Import GO terms and orthologs to Calc
# 5) Use the *.bed.clean files to pair-up hits with genes and GO terms







#!/bin/bash
# IEKNGHER1ujm

# ==============================================================================
# Genotype haploids with snps-diploid.sh to recover heterozygous SNP calls
# Het SNPs that are close (~2kb) are potential CNVs and may be supported by
# high depth of coverage. Purpose of this script is to investigate this.
# ==============================================================================

# Each sample needs processing individually because:
# 1) Depth of coverage will be important and not all samples are "equal" here
# 2) Imputed genotypes should not be used

# Plan is as follows:
# 1) Genotype haploids using the snps-diploid.sh script
# 2) Extract SNPs with heterozygous genotypes
# 3) Analyse distribution and frequency of SNPs
# 4) Analyse depth of coverage with respect to sample average
# Load modules

module load bioinfo/bcftools

# Paths to main folders
DUMP=/home/dwragg/work/Analysis
PIPE=/save/seqapipop/Scripts
PERL=${PIPE}/Perl
SEQA=/save/seqapipop/Data/Apis-mellifera
AMEL=/work/dwragg/Apis-mellifera
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"



# Script is now stored as scrib-hetSNPs.sh
qsub -q workq -l mem=8G -l h_vmem=48G -pe parallel_smp 1 \
    -o ${DUMP}/hetSNPs -e ${DUMP}/hetSNPs ${PIPE}/scrib-hetSNPs.sh


# Plink
${PLINK}/plink --vcf ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf \
  --keep-allele-order \
  --a2-allele ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf 4 3 '#' \
  --allow-no-sex \
  --allow-extra-chr \
  --mind 1 \
  --geno 0.9 \
  --maf 0.1 \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/hetSNPs/JFMOTH-hetSNPs \
  --freq \
  --missing \
  --make-bed

bgzip -f ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf
tabix -fp vcf ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz
gunzip -c ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz > ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf

# Record for each sample a distant snplist, necessary because indications of
# a CNV within A SAMPLE will be 2+ het SNPs within the specified distance
POP="JFM"
SAMPLES=($(cat ${DUMP}/vcfs/lists/${POP}_samples_fixed.list))
for SAMPLE in ${SAMPLES[@]}
do
  echo ${SAMPLE} > ${DUMP}/hetSNPs/sample

  # Identify SNPs outside 2kb
  ${PLINK}/plink --bfile ${DUMP}/hetSNPs/JFMOTH-hetSNPs \
    --keep-fam ${DUMP}/hetSNPs/sample \
    --keep-allele-order \
    --a2-allele ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf 4 3 '#' \
    --allow-no-sex \
    --allow-extra-chr \
    --mind 1 \
    --geno 0 \
    --maf 0.1 \
    --chr-set 16 \
    --chr 1-16 \
    --set-missing-snp-ids @:#[Amel4-5] \
    --out ${DUMP}/hetSNPs/snplists/${SAMPLE}-distant \
    --bp-space 2000 \
    --write-snplist
  rm ${DUMP}/hetSNPs/snplists/${SAMPLE}-distant.log
  rm ${DUMP}/hetSNPs/snplists/${SAMPLE}-distant.nosex

  # Filter out the distant SNPs and generate a snplist
  ${PLINK}/plink --bfile ${DUMP}/hetSNPs/JFMOTH-hetSNPs \
    --keep-fam ${DUMP}/hetSNPs/sample \
    --exclude ${DUMP}/hetSNPs/snplists/${SAMPLE}-distant.snplist \
    --keep-allele-order \
    --a2-allele ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf 4 3 '#' \
    --allow-no-sex \
    --allow-extra-chr \
    --mind 1 \
    --geno 0 \
    --maf 0.1 \
    --chr-set 16 \
    --chr 1-16 \
    --set-missing-snp-ids @:#[Amel4-5] \
    --out ${DUMP}/hetSNPs/vcfs/${SAMPLE} \
    --recode vcf-iid
  rm ${DUMP}/hetSNPs/vcfs/${SAMPLE}.log
  rm ${DUMP}/hetSNPs/vcfs/${SAMPLE}.nosex

  # Extract SNPs using bcftools to retain DP
  bcftools view -R ${DUMP}/hetSNPs/vcfs/${SAMPLE}.vcf -s ${SAMPLE} \
    -O v -o ${DUMP}/hetSNPs/vcfs/${SAMPLE}.vcf \
    ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz

  # Generate bed file from VCF
  awk '{OFS="\t"; if(!/^#/){ split($10, a, ":"); print $1,$2,$2,a[2]}}' \
    ${DUMP}/hetSNPs/vcfs/${SAMPLE}.vcf \
    > ${DUMP}/hetSNPs/beds/${SAMPLE}.bed

done



POP="JFM"
SAMPLES=($(cat ${DUMP}/vcfs/lists/${POP}_samples_fixed.list))
for SAMPLE in ${SAMPLES[@]}
do
  bedtools merge -d 2000  -i ${DUMP}/hetSNPs/beds/${SAMPLE}.bed \
    -c 4,4 -o count,mean \
    > ${DUMP}/hetSNPs/beds/${SAMPLE}.merge.bed
done













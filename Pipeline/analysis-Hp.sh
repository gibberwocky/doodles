#!/bin/bash
# IEKNGHER1ujm

# Set-up
module load bioinfo/bcftools
VCFT=/usr/local/bioinfo/src/vcftools/current/perl
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
SITES=${DUMP}/vcfs/HARP-JFM-OTH-07112014



# Use BED file from SITES
OUT=${DUMP}/selection/hp
FILE=${DUMP}/selection/S2_hapJFMOTH-plink
PREFIX="S2-07112014"
POP="JFM"

bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3' \
-S ${DUMP}/vcfs/lists/${POP}_samples.list \
  -o ${OUT}/${PREFIX}-${POP}.vcf \
  ${SITES}/S2_hapJFMOTH_beagle.vcf.gz

# Change from accession to chromosome naming
${PIPE}/accession-to-chr.sh -i ${OUT}/${PREFIX}-${POP}.vcf
sed -i '/^#CHROM/s/_/-/g' ${OUT}/${PREFIX}-${POP}.vcf


# Perform test for signature of selection using pooled heterozygosity in
# sliding windows along a chromosome.

# =============================================================================
# Perl (quickest method of the 3 listed here)
#	- pooled heterozygosity approach (as per Rubin et al in chicken)
# =============================================================================

SCRIPTS="/save/dwragg/Scripts/Perl"

# Fill AN if missing:
perl ${VCFT}/fill-an-ac ${OUT}/${PREFIX}-${POP}.vcf > ${OUT}/${PREFIX}-${POP}-AN.vcf
mv -f ${OUT}/${PREFIX}-${POP}-AN.vcf ${OUT}/${PREFIX}-${POP}.vcf


# Windowed Hp analysis
perl ${SCRIPTS}/VCF-Hp.pl ${OUT}/${PREFIX}-${POP}.vcf 16 5000 1000


# Single Hp value genome-wide (for "heterozygosity")
perl ${SCRIPTS}/VCF-Hp-global.pl ${OUT}/${PREFIX}-${POP}.vcf 

# Check heterozygosity in AMCO lineages
POP="O"
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3' \
  -o ${OUT}/${POP}.vcf \
  ${SITES}/${POP}_beagle.vcf.gz
# Fill AN if missing:
perl ${VCFT}/fill-an-ac ${OUT}/${POP}.vcf > ${OUT}/${POP}-AN.vcf
# Calculate Hp
perl ${SCRIPTS}/VCF-Hp-global.pl ${OUT}/${POP}-AN.vcf


# =============================================================================
# VCFTools + R
#	- pooled heterozygosity approach (as per Rubin et al in chicken)
# =============================================================================

# Generate allele frequency counts
vcftools --vcf ${OUT}/${PREFIX}-${POP}.vcf \
  --out ${OUT}/${PREFIX}-${POP} --counts

qsub -q unlimitq -l mem=4G -l h_vmem=32G -pe parallel_smp 1 \
  -o ${DUMP}/selection/hp/logs \
  -e ${DUMP}/selection/hp/logs \
  ${PIPE}/Hp-VCFTools-counts.sh \
    -d ${OUT} \
    -i ${PREFIX}-${POP}.frq.count \
    -w 5000 \
    -x 1000




# =============================================================================
# R - (alternative method)
#	- pooled heterozygosity approach (as per Rubin et al in chicken)
# =============================================================================

# Loop through each autosome and create files for R script
mkdir -p ${DUMP}/selection/hp/logs
mkdir -p ${DUMP}/selection/hp/${PREFIX}-${POP}
cd ${DUMP}/selection/hp
for i in {1..16}
do
  # Create input file containing only chr of interest
  ${PLINK}/plink --vcf ${OUT}/${PREFIX}-${POP}.vcf \
    --allow-no-sex \
    --allow-extra-chr \
    --chr-set 16 \
    --chr ${i} \
    --set-missing-snp-ids @:#[Amel4-5] \
    --recode 12 \
    --out ${DUMP}/selection/hp/${PREFIX}-${POP}/${i}
  # Remove genetic distance from map file
  cut -f1,2,4 ${DUMP}/selection/hp/${PREFIX}-${POP}/${i}.map > \
    ${DUMP}/selection/hp/${PREFIX}-${POP}/${i}-R.map
  rm ${DUMP}/selection/hp/${PREFIX}-${POP}/${i}.log
  rm ${DUMP}/selection/hp/${PREFIX}-${POP}/${i}.nosex
  rm ${DUMP}/selection/hp/${PREFIX}-${POP}/${i}.map
done

# Loop through each autosome and run pooled het script
for i in {1..16}
do
  # Run pooled heterozygosity script
  qsub -q unlimitq -l mem=4G -l h_vmem=32G -pe parallel_smp 1 \
      -o ${DUMP}/selection/hp/logs \
      -e ${DUMP}/selection/hp/logs \
      ${PIPE}/selection.sh \
        -i ${DUMP}/selection/hp/${PREFIX}-${POP}/ \
        -p ${DUMP}/selection/hp/${PREFIX}-${POP}/${i}.ped \
        -m ${DUMP}/selection/hp/${PREFIX}-${POP}/${i}-R.map \
        -h ${DUMP}/selection/hp/S2_hapJFMOTH-plink.pheno \
        -c ${i} \
        -w 5000 \
        -o ${DUMP}/selection/hp/${PREFIX}-${POP} \
	-x 1000
done



# Remove obsolete files, and merge results
# Only run once above iterations completed!
# Make sure filename correct for *.txt because R sometimes using
# scientific notation for large values (window sizes)
POP="OTH"
for i in {1..16}
do
  # Remvoe obsolte files
#  rm ${DUMP}/selection/${pop}/${i}.ped
#  rm ${DUMP}/selection/${pop}/${i}-R.map
#  rm ${DUMP}/selection/${pop}/gen${i}

  # Merge results
  cd ${DUMP}/selection/hp/${PREFIX}-${POP}
  # Ensure -name is correct
  find . -name "*.txt" | xargs -n 1 tail -n +2 | \
    cut -f 2- > ${DUMP}/selection/hp/${PREFIX}-${POP}/Hp.body
  # ensure file name is correct
  head --l 1 ${DUMP}/selection/hp/${PREFIX}-${POP}/*chr1-*.txt > \
    ${DUMP}/selection/hp/${PREFIX}-${POP}/Hp.head
  cat ${DUMP}/selection/hp/${PREFIX}-${POP}/Hp.head \
    ${DUMP}/selection/hp/${PREFIX}-${POP}/Hp.body > \
    ${DUMP}/selection/hp/${PREFIX}-${POP}.out
#  rm ${DUMP}/selection/hp/${POP}/Hp.head
#  rm ${DUMP}/selection/hp/${POP}/Hp.body
#  rm ${DUMP}/selection/hp/${POP}/*.txt
done
















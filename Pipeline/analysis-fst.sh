#!/bin/bash
# IEKNGHER1ujm


# Set-up
module load bioinfo/bcftools
VCFT="/usr/local/bioinfo/src/vcftools/vcftools_0.1.12a/bin"
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
SITES=${DUMP}/vcfs/HARP-JFM-OTH-07112014
OUTFILE="S2-07112014_hapJFMOTH"

# Reduce to chromosomes 1-16
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3' \
  -o ${DUMP}/selection/vcft-fst/${OUTFILE}.vcf \
  ${SITES}/S2_hapJFMOTH_beagle.vcf.gz



# Change from accession to chromosome naming
${PIPE}/accession-to-chr.sh -i ${DUMP}/selection/vcft-fst/${OUTFILE}.vcf 
sed -i '/^#CHROM/s/_/-/g' ${DUMP}/selection/vcft-fst/${OUTFILE}.vcf 


cp ${DUMP}/vcfs/lists/JFM_samples.list ${DUMP}/selection/vcft-fst/JFM_samples.list
sed -i 's/_/-/g' ${DUMP}/selection/vcft-fst/JFM_samples.list
cp ${DUMP}/vcfs/lists/OTH_samples.list ${DUMP}/selection/vcft-fst/OTH_samples.list
sed -i 's/_/-/g' ${DUMP}/selection/vcft-fst/OTH_samples.list

# =============================================================================
# VCFtools
#	- windowed Fst between populations (Weir & Cockerham 1984)
# =============================================================================

# Calculate windowed Fst (Weir & Cockerham 1984)
${VCFT}/vcftools --vcf ${DUMP}/selection/vcft-fst/${OUTFILE}.vcf \
  --weir-fst-pop ${DUMP}/selection/vcft-fst/JFM_samples.list \
  --weir-fst-pop ${DUMP}/selection/vcft-fst/OTH_samples.list \
  --fst-window-size 5000 \
  --fst-window-step 1000 \
  --out ${DUMP}/selection/vcft-fst/${OUTFILE}

# 		window-size	window-step	u(SNPs)
#	T4	5000		1000		83
#	T4-2	10000		2000		164
#	T4-3	40000		10000		627

# * too few SNPs

# Output summary stats
awk 'NR == 1 { sum=0 } { sum+=$4;} END {printf "Average: %f\n", sum/NR}' \
  ${DUMP}/selection/vcft-fst/${OUTFILE}*.fst





# =============================================================================
# R - Selection 
#	- windowed Fst between populations
# =============================================================================

# Performs test for signature of selection by calculating Fst and Josts' D in
# sliding windows along a chromosome.

# Rec-cap on parameters
DUMP="/home/dwragg/work/Analysis"
PIPE="/home/dwragg/work/Pipeline"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
OUT=${DUMP}/selection/R-fst
FILE=${OUT}/S2_hapJFMOTH-plink
pop="fst"


# Make BED file
${PLINK}/plink --vcf ${DUMP}/selection/vcft-fst/S2_hapJFMOTH.vcf \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${FILE} \
  --make-bed

# FAM to PHENO (id	pop	sex	 phene)
# manually edit pop field
#cut -d' ' -f1,3,4 ${FILE}.fam > ${FILE}.pheno


# Loop through each autosome and create files for R script
mkdir -p ${OUT}/${pop}/logs
cd ${OUT}/${pop}
for i in {1..16}
do
  # Create input file containing only chr of interest
  ${PLINK}/plink --bfile ${FILE} \
    --allow-no-sex \
    --allow-extra-chr \
    --chr-set 16 \
    --chr ${i} \
    --set-missing-snp-ids @:#[Amel4-5] \
    --recode 12 \
    --out ${OUT}/${pop}/${i}
  # Remove genetic distance from map file
  cut -f1,2,4 ${OUT}/${pop}/${i}.map > ${OUT}/${pop}/${i}-R.map
  rm ${OUT}/${pop}/${i}.log
  rm ${OUT}/${pop}/${i}.nosex
  rm ${OUT}/${pop}/${i}.map
done

# Loop through each autosome and run winFst script
for i in {1..16}
do
  # Run pooled heterozygosity script
  qsub -q unlimitq -l mem=4G -l h_vmem=32G -pe parallel_smp 1 \
      -o ${OUT}/${pop}/logs \
      -e ${OUT}/${pop}/logs \
      ${PIPE}/winfst.sh \
        -i ${DUMP}/vcfs \
        -p ${OUT}/${pop}/${i}.ped \
        -m ${OUT}/${pop}/${i}-R.map \
        -h ${FILE}.pheno \
        -c ${i} \
        -w 5000 \
        -o ${OUT}/${pop} \
	-x 1000
done

# ============================================================================
# <-         upt to here R script needs to be checked to ensure QC is correct
# ============================================================================

# Remove obsolete files, and merge results
# Only run once above iterations completed!
#rename -5000 -5k ${OUT}/${pop}/*
for i in {1..16}
do
  # Remvoe obsolte files
#  rm ${OUT}/${pop}/gen${i}
  # Merge results
  cd ${OUT}/${pop}/
  find . -name "*.txt" | xargs -n 1 tail -n +2 | cut -f 2- > ${OUT}/${pop}/Fst.body
  head --l 1 ${OUT}/${pop}/*chr1-*.txt > ${OUT}/${pop}/Fst.head
  cat ${OUT}/${pop}/Fst.head ${OUT}/${pop}/Fst.body > ${OUT}/R-fst.out
#  rm ${OUT}/${pop}/Fst.head
#  rm ${OUT}/${pop}/Fst.body
done





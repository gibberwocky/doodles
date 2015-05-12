#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

# Paths
DUMP="/home/dwragg/work/Analysis"
SITES=${DUMP}/vcfs/HARP-JFM-OTH-07112014
FILE=Admixture-dip-prune
PIPE="/save/seqapipop/Scripts"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
VCFT="/usr/local/bioinfo/src/vcftools/vcftools_0.1.12a/bin"

# Write snps
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE} \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --recode vcf \
  --out ${DUMP}/vcfs/plink/${FILE}-tmp

# Translate chromosome to NCBI accesions
qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/vcfchr-to-accession2.sh -i ${DUMP}/vcfs/plink/${FILE}-tmp.vcf

# Cut chromosome and position
cut -f1,2 ${DUMP}/vcfs/plink/${FILE}-tmp.vcf > ${DUMP}/vcfs/plink/${FILE}-tmp.bed
sed -i '/^#/ d' ${DUMP}/vcfs/plink/${FILE}-tmp.bed
sort -k1,1 -k2n ${DUMP}/vcfs/plink/${FILE}-tmp.bed > ${DUMP}/vcfs/plink/${FILE}-tmp2.bed

# Output site details to Admixture folder
sites=${DUMP}/sites/HARP-JFM-OTH.vcf.gz
out=${DUMP}/vcfs/plink/${FILE}.vcf
cd ${DUMP}/vcfs/plink
${VCFT}/vcftools --gzvcf ${sites} \
  --positions ${DUMP}/vcfs/plink/${FILE}-tmp2.bed \
  --recode
mv ${DUMP}/vcfs/plink/out.recode.vcf ${DUMP}/vcfs/plink/${FILE}.vcf
rm ${DUMP}/vcfs/plink/out.log
rm ${DUMP}/vcfs/plink/${FILE}*tmp*
bgzip -f ${DUMP}/vcfs/plink/${FILE}.vcf
tabix -fp vcf ${DUMP}/vcfs/plink/${FILE}.vcf.gz









# ==============================================================================
# Quick test of newly sequenced data for ADMIXTURE ancestry assignment
# ==============================================================================

# Output path and master site file
sites=${DUMP}/vcfs/plink/${FILE}.vcf.gz
OUT=${DUMP}/vcfs/plink
POP="O"

# Run each dataset (A, M, C, O, etc) on cluster (check ploidy -p)
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/imp-geno.sh -o ${OUT} -s ${sites} -p 2 \
      -v ${DUMP}/vcfs/lists/${POP}_vcfs.list \
      -b ${DUMP}/vcfs/lists/${POP}_bams.list \
      -n ${POP}-33k

# Clean up unnecessary files
rm ${OUT}/*subhet*



# Set-up
module load bioinfo/bcftools
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
SITES=${DUMP}/vcfs/plink
FILEDIR=${DUMP}/vcfs/VCF-SITES

# Root file name for Plink/Admixture input/output
FILE=Admixture-dip

# Merge data (diploids for admixture)
bcftools merge -m snps -O z \
  ${SITES}/A-Harpur_beagle.vcf.gz \
  ${SITES}/M-Harpur_beagle.vcf.gz \
  ${SITES}/C-Harpur_beagle.vcf.gz \
  ${SITES}/O-Harpur_beagle.vcf.gz \
  ${SITES}/JFM_diploid.vcf.gz \
  ${SITES}/OTH_diploid.vcf.gz \
  ${SITES}/BR_diploid.vcf.gz \
  ${SITES}/OUE_diploid.vcf.gz \
  ${SITES}/COR_diploid.vcf.gz \
  -o ${SITES}/test-Harpur_beagle.vcf.gz
tabix -fp vcf ${SITES}/test-Harpur_beagle.vcf.gz

# Root file name
INFILE=test-Harpur_beagle.vcf

# Reduce to chromosomes 1-16 + Mt
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3','NC_001566.1' \
  -o ${FILEDIR}/${FILE}.vcf \
  ${SITES}/${INFILE}.gz

# Replace underscores in sample IDs
sed -i '/^#CHROM/s/_/-/g' ${FILEDIR}/${FILE}.vcf
# Convert chromosome accessions to numbers
${PIPE}/accession-to-chr.sh -i ${FILEDIR}/${FILE}.vcf

# Make BED file
${PLINK}/plink --vcf ${FILEDIR}/${FILE}.vcf \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --make-bed \
  --out ${DUMP}/vcfs/plink/${FILE}-prune

# Run ADMIXTURE
k=1 # Min K
K=8 # Max K
N=1 # Number of bootstraps
# Check the -i file is correct (pruned or not) before running
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/admixture.sh -i ${DUMP}/vcfs/plink/${FILE}-prune \
      -k ${k} -K ${K} -n ${N} -o ${DUMP}/vcfs/admixture -s ${FILE}-prune
# Identify optimal K
grep -h CV ${DUMP}/vcfs/admixture/${FILE}-*.out
# Copy fam file to same folder incase needed for plotting
cp ${DUMP}/vcfs/plink/${FILE}-prune.fam ${DUMP}/vcfs/admixture









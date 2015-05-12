#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

# Paths
PIPE="/save/seqapipop/Scripts"
DUMP="/home/dwragg/work/Analysis"
SITES="${DUMP}/vcfs/HARP-JFM-OTH-07112014"
FILE="S2_AMCO.vcf"
FILEDIR="${DUMP}/vcfs/VCF-SITES"
FILE1="Harpur"
VCFT="/usr/local/bioinfo/src/vcftools/vcftools_0.1.12a/bin"


# ==============================================================================
# 1. Make VCF file from Harpur data for filtering with Plink
# ==============================================================================

# Merge data (reference only for admixture)
bcftools merge -m snps -O z \
  ${SITES}/A_beagle.vcf.gz \
  ${SITES}/M_beagle.vcf.gz \
  ${SITES}/C_beagle.vcf.gz \
  ${SITES}/O_beagle.vcf.gz \
  -o ${SITES}/${FILE}.gz
tabix -fp vcf ${SITES}/${FILE}.gz

# Reduce to chromosomes 1-16 + Mt
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3','NC_001566.1' \
  -o ${FILEDIR}/${FILE1}.vcf \
  ${SITES}/${FILE}.gz

# Replace underscores in sample IDs
sed -i '/^#CHROM/s/_/-/g' ${FILEDIR}/${FILE1}.vcf
# Convert chromosome accessions to numbers
${PIPE}/accession-to-chr.sh -i ${FILEDIR}/${FILE1}.vcf


# ==============================================================================
# 2. Filter Harpur VCF file with Plink
# ==============================================================================

# Make BED file
${PLINK}/plink --vcf ${FILEDIR}/${FILE1}.vcf \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/plink/${FILE1} \
  --make-bed

# Before running ADMIXTURE first remove SNPs in high LD
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE1} \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/plink/${FILE1} \
  --indep-pairwise 50 10 0.1
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE1} \
  --allow-no-sex \
  --allow-extra-chr \
  --geno 0.1 \
  --maf 0.01 \
  --bp-space 5000 \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --extract ${DUMP}/vcfs/plink/${FILE1}.prune.in \
  --out ${DUMP}/vcfs/plink/${FILE1}-prune \
  --make-bed


# ==============================================================================
# 3. Recover filtered SNPs to create reference panel for genotyping
# ==============================================================================

# Write snps
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE1}-prune \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --recode vcf \
  --out ${FILEDIR}/${FILE1}-prune-SNPs

# Translate chromosome to NCBI accesions
qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/vcfchr-to-accession2.sh -i ${FILEDIR}/${FILE1}-prune-SNPs.vcf

# Cut chromosome and position
cut -f1,2 ${FILEDIR}/${FILE1}-prune-SNPs.vcf > ${FILEDIR}/${FILE1}-prune-SNPs.tmp
sed -i '/^#/ d' ${FILEDIR}/${FILE1}-prune-SNPs.tmp
sort -k1,1 -k2n ${FILEDIR}/${FILE1}-prune-SNPs.tmp > ${FILEDIR}/${FILE1}-prune-SNPs.bed
rm ${FILEDIR}/${FILE1}-prune-SNPs.tmp

# Output site details to Admixture folder
out=${DUMP}/vcfs/plink/${FILE}.vcf
cd ${FILEDIR}
${VCFT}/vcftools --gzvcf ${SITES}/${FILE}.gz \
  --positions ${FILEDIR}/${FILE1}-prune-SNPs.bed \
  --recode
rm ${FILEDIR}/out.log
rm ${FILEDIR}/${FILE1}* 
mv ${FILEDIR}/out.recode.vcf ${FILEDIR}/${FILE1}-prune.vcf
bgzip -f ${FILEDIR}/${FILE1}-prune.vcf
tabix -fp vcf ${FILEDIR}/${FILE1}-prune.vcf.gz


# ==============================================================================
# 4. Genotype all samples for SNPs in the reference panel
# ==============================================================================

# Output path and master site file
sites=${FILEDIR}/${FILE1}-prune.vcf.gz
OUT=${DUMP}/vcfs/plink
POP="COR"

# Run each dataset (A, M, C, O, etc) on cluster (check ploidy -p)
qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/imp-geno.sh -o ${OUT} -s ${sites} -p 1 \
      -v ${DUMP}/vcfs/lists/${POP}_vcfs.list \
      -b ${DUMP}/vcfs/lists/${POP}_bams.list \
      -n ${POP}-Harpur

# Clean up unnecessary files once above completed for all samples
rm ${OUT}/*subhet*


# ==============================================================================
# 5. Create diploids from haploids (requires even number of samples as input!)
# ==============================================================================

# Identify samples for diploidy
SITES=${DUMP}/vcfs/plink
POP="COR"
SUFFIX="-Harpur_beagle"
readarray HAPS < ${DUMP}/vcfs/lists/${POP}_samples.list

# Create file to record sample IDs for merging
rm -f ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
touch ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
# Loop through all pairs
x=0
for ((i=0; i<${#HAPS[@]}; i=i+2))
do
  # Write sample ID to file
  printf "${SITES}/%s.vcf.gz\n" "${POP}-DIP${x}" \
    >> ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
  # Create directory for logs
  mkdir -p ${DUMP}/logs/${POP}-DIP${x}
  # Create diploid from samples i and i+1
  qsub -q workq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/${POP}-DIP${x} \
    -e ${DUMP}/logs/${POP}-DIP${x} \
    ${PIPE}/diploid-drone.sh \
      -v ${SITES}/${POP}${SUFFIX}.vcf.gz \
      -a ${HAPS[${i}]} \
      -b ${HAPS[${i}+1]} \
      -s ${POP}-DIP${x} \
      -o ${SITES}
  x=$((x+1))
done

# Merge diploids into single VCF
bcftools merge -m snps -O z \
  -l ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list \
  -o ${SITES}/${POP}_diploid.vcf.gz
tabix -p vcf ${SITES}/${POP}_diploid.vcf.gz
# Remove unnecessary files
readarray FILES < ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
rm ${FILES[@]}
tmp=$(printf "%s.tbi " ${FILES[@]})
rm ${tmp}


# ==============================================================================
# 6. Run ADIMXTURE
# ==============================================================================

# Variables
FILEDIR=${DUMP}/vcfs/VCF-SITES
INFILE=S3_AMCO-dipALL.vcf
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
  ${SITES}/BER_diploid.vcf.gz \
  -o ${SITES}/${INFILE}.gz
tabix -fp vcf ${SITES}/${INFILE}.gz


# Reduce to chromosomes 1-16 + Mt
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3','NC_001566.1' \
  -o ${SITES}/${FILE}.vcf \
  ${SITES}/${INFILE}.gz

# Replace underscores in sample IDs
sed -i '/^#CHROM/s/_/-/g' ${SITES}/${FILE}.vcf
# Convert chromosome accessions to numbers
${PIPE}/accession-to-chr.sh -i ${SITES}/${FILE}.vcf

# Make BED file
${PLINK}/plink --vcf ${SITES}/${FILE}.vcf \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/plink/${FILE} \
  --make-bed

k=1 # Min K
K=8 # Max K
N=1 # Number of bootstraps
# Check the -i file is correct (pruned or not) before running
qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/admixture.sh -i ${DUMP}/vcfs/plink/${FILE} \
      -k ${k} -K ${K} -n ${N} -o ${DUMP}/vcfs/admixture -s ${FILE}

# Identify optimal K
grep -h CV ${DUMP}/vcfs/admixture/${FILE}-*.out
# Copy fam file to same folder incase needed for plotting
cp ${DUMP}/vcfs/plink/${FILE}.fam ${DUMP}/vcfs/admixture




# =============================================================================
# 7. Phylogenetic analysis
# =============================================================================

# This assumes that the Plink file used in the ADMIXTURE analysis is intacted
# and in the ${DUMP}/vcfs/plink folder, and that the ADMIXTURE results have
# been transferred to their own folder in ${DUMP}/vcfs/admixture, including
# the plink fam file.

INDIR="admixture/S3-All"
FILE="Admixture-dip"

# FAM to PHENO (id	sex	 phene)
cut -d' ' -f1,3,4 ${DUMP}/vcfs/${INDIR}/${FILE}.fam \
  > ${DUMP}/vcfs/${INDIR}/${FILE}-R.tmp
header=('id' 'sex' 'pop')
echo ${header[@]} | cat - ${DUMP}/vcfs/${INDIR}/${FILE}-R.tmp \
  > ${DUMP}/vcfs/${INDIR}/${FILE}-R.pheno
rm ${DUMP}/vcfs/${INDIR}/${FILE}-R.tmp

# Format Plink for GenABEL import
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE} \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --recode 12 \
  --out ${DUMP}/vcfs/${INDIR}/${FILE}-R

# Remove genetic distance from map file
cut -f1,2,4 ${DUMP}/vcfs/${INDIR}/${FILE}-R.map \
  > ${DUMP}/vcfs/${INDIR}/${FILE}-R.tmp
mv ${DUMP}/vcfs/${INDIR}/${FILE}-R.tmp ${DUMP}/vcfs/${INDIR}/${FILE}-R.map

# Proceed in R




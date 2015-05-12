#!/bin/bash
# IEKNGHER1ujm

# Set-up
module load bioinfo/bcftools
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
#SITES=${DUMP}/vcfs/AMEL-4M
SITES=${DUMP}/vcfs/HARP-JFM-OTH-07112014
FILEDIR=${DUMP}/vcfs/VCF-SITES

# Root file name for Plink/Admixture input/output
FILE=Admixture-dip



# *Consider using just a set of SNPs identified in Harpur's data


# =============================================================================
# Pre-process VCF for Plink/Admixture
# =============================================================================

# Merge data (reference only for admixture)
#bcftools merge -m snps -O z \
#  ${SITES}/A_beagle.vcf.gz \
#  ${SITES}/M_beagle.vcf.gz \
#  ${SITES}/C_beagle.vcf.gz \
#  ${SITES}/O_beagle.vcf.gz \
#  -o ${SITES}/S1_AMCO.vcf.gz
#tabix -fp vcf ${SITES}/S2_AMCO.vcf.gz

# Merge data (diploids for admixture)
bcftools merge -m snps -O z \
  ${SITES}/A_beagle.vcf.gz \
  ${SITES}/M_beagle.vcf.gz \
  ${SITES}/C_beagle.vcf.gz \
  ${SITES}/O_beagle.vcf.gz \
  ${SITES}/JFM_diploid.vcf.gz \
  ${SITES}/OTH_diploid.vcf.gz \
  -o ${SITES}/S2_AMCO-dipJFMOTH_beagle.vcf.gz
tabix -fp vcf ${SITES}/S2_AMCO-dipJFMOTH_beagle.vcf.gz

# Merge data (diploids for admixture with Wallberg)
#bcftools merge -m snps -O z \
#  ${SITES}/A_HW_beagle.vcf.gz \
#  ${SITES}/M_HW_beagle.vcf.gz \
#  ${SITES}/C_HW_beagle.vcf.gz \
#  ${SITES}/O_HW_beagle.vcf.gz \
#  ${SITES}/JFM_diploid.vcf.gz \
#  ${SITES}/OTH_diploid.vcf.gz \
#  -o ${SITES}/S2_AMCOHW-dipJFMOTH_beagle.vcf.gz
#tabix -fp vcf ${SITES}/S2_AMCOHW-dipJFMOTH_beagle.vcf.gz

# Merge data (haploids for admixture)
#bcftools merge -m snps -O z \
#  ${SITES}/A_beagle.vcf.gz \
#  ${SITES}/M_beagle.vcf.gz \
#  ${SITES}/C_beagle.vcf.gz \
#  ${SITES}/O_beagle.vcf.gz \
#  ${SITES}/JFM_beagle.vcf.gz \
#  ${SITES}/OTH_beagle.vcf.gz \
#  -o ${SITES}/S2_AMCO-hapJFMOTH.vcf.gz
#tabix -fp vcf ${SITES}/S2_AMCO-hapJFMOTH.vcf.gz

# Root file name
INFILE=S2_AMCO-dipJFMOTH_beagle.vcf

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


# =============================================================================
# PLINK
# =============================================================================

# Make BED file
${PLINK}/plink --vcf ${FILEDIR}/${FILE}.vcf \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/plink/${FILE} \
  --make-bed


# Before running ADMIXTURE first remove SNPs in high LD
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE} \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/plink/${FILE} \
  --indep-pairwise 50 10 0.1
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE} \
  --allow-no-sex \
  --allow-extra-chr \
  --geno 0.1 \
  --maf 0.01 \
  --bp-space 5000 \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --extract ${DUMP}/vcfs/plink/${FILE}.prune.in \
  --out ${DUMP}/vcfs/plink/${FILE}-prune \
  --make-bed
# bp-space retains only 1 SNP from a pair closer than the given bp
# (Harpur et al used 25K random SNPs at least 5 kb apart)
# (Wallberg used all 8M)
# (Here I've used LD cut-off as per Admixture manual, that is all)

# Optional prune on chromosome also
CHR=16
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE}-prune \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr ${CHR} \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/plink/${FILE}-${CHR} \
  --make-bed


# Check fraction missing per sample (update to check before or after pruning)
# To look at an individual chromosome change the --chr value
#${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE}-prune \
#  --allow-no-sex \
#  --allow-extra-chr \
#  --chr-set 16 \
#  --chr 1-16 \
#  --set-missing-snp-ids @:#[Amel4-5] \
#  --out ${DUMP}/vcfs/plink/${FILE}-prune \
#  --missing
# Output average of col 5 to save importing to Calc
#awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' ${DUMP}/vcfs/plink/${FILE}-prune.imiss


# =============================================================================
# ADMIXTURE
# =============================================================================

# Runs ADMIXTURE for several K values (1 to K)
# The lowest cross validation (CV) error rate is the optimal K
# Input format is plink
# ADMIXTURE outputs file for:
#  Q (ancestry fractions)
#  P (allele frequencies of inferred ancestral populations)
k=1 # Min K
K=8 # Max K
N=1 # Number of bootstraps
CHR=2 # if not chromosome-specific then change CHR to prune
PREFIX=S2-CHR-${CHR}
# Check the -i file is correct (pruned or not) before running
mkdir -p ${DUMP}/vcfs/admixture/${PREFIX}
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/admixture.sh -i ${DUMP}/vcfs/plink/${FILE}-${CHR} \
      -k ${k} -K ${K} -n ${N} -o ${DUMP}/vcfs/admixture/${PREFIX}/ -s ${FILE}-${CHR}
# Copy fam file to same folder incase needed for plotting
cp ${DUMP}/vcfs/plink/${FILE}-${CHR}.fam ${DUMP}/vcfs/admixture${PREFIX}
# Identify optimal K
grep -h CV ${DUMP}/vcfs/admixture/${PREFIX}/${FILE}-*.out





# =============================================================================
# Plink to R - for phylogenetic analysis
# =============================================================================

# This assumes that the Plink file used in the ADMIXTURE analysis is intacted
# and in the ${DUMP}/vcfs/plink folder, and that the ADMIXTURE results have
# been transferred to their own folder in ${DUMP}/vcfs/admixture, including
# the plink fam file.

INDIR="admixture/S2-Harpur-OUECOR"
FILE="Admixture-dip"

# FAM to PHENO (id	sex	 phene)
cut -d' ' -f1,3,4 ${DUMP}/vcfs/${INDIR}/${FILE}-prune.fam \
  > ${DUMP}/vcfs/${INDIR}/${FILE}-prune-R.tmp
header=('id' 'sex' 'pop')
echo ${header[@]} | cat - ${DUMP}/vcfs/${INDIR}/${FILE}-prune-R.tmp \
  > ${DUMP}/vcfs/${INDIR}/${FILE}-prune-R.pheno
rm ${DUMP}/vcfs/${INDIR}/${FILE}-prune-R.tmp

# Format Plink for GenABEL import
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${FILE}-prune \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --recode 12 \
  --out ${DUMP}/vcfs/${INDIR}/${FILE}-prune-R

# Remove genetic distance from map file
cut -f1,2,4 ${DUMP}/vcfs/${INDIR}/${FILE}-prune-R.map \
  > ${DUMP}/vcfs/${INDIR}/${FILE}-prune-R.tmp
mv ${DUMP}/vcfs/${INDIR}/${FILE}-prune-R.tmp ${DUMP}/vcfs/${INDIR}/${FILE}-prune-R.map

# Proceed in R





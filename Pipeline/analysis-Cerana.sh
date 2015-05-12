#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

# Paths to main folder
ACER=/work/dwragg/Apis-cerana
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
OUT=${DUMP}/vcfs/HARP-JFM-OTH-15012015
sites=${DUMP}/sites/Apis-cerana

#SEQA=/save/seqapipop/Data/Apis-cerana/Harpur
#mkdir ${ACER}/vcfs
#mkdir ${ACER}/bams
#ln -fs ${SEQA}/*/vcfs/*_clean.vcf.gz ${ACER}/vcfs/
#ln -fs ${SEQA}/*/*_bootstrap* ${ACER}/bams/
# tabix -fp vcf ${ACER}/vcfs/SRR957079_clean.vcf.gz

# Subset Apis-cerana SNPs to those with good DP quality
bcftools view --types snps \
  --include 'DP>=9 & DP<(AVG(DP)*3)' \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3' \
  -M 2 -O v -o ${OUT}/Apis-cerana.vcf \
  ${ACER}/vcfs/SRR957079_clean.vcf.gz

# Change from accession to chromosome naming and fix _ in sample names
${PIPE}/accession-to-chr.sh -i ${OUT}/Apis-cerana.vcf
sed -i '/^#CHROM/s/_/-/g' ${OUT}/Apis-cerana.vcf

# Feed through Plink to populate SNP IDs
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
${PLINK}/plink --vcf ${OUT}/Apis-cerana.vcf \
  --keep-allele-order \
  --a2-allele ${OUT}/Apis-cerana.vcf 4 3 '#' \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${OUT}/Apis-cerana_QC \
  --recode vcf-iid
mv ${OUT}/Apis-cerana_QC.vcf ${OUT}/Apis-cerana.vcf
rm ${OUT}/Apis-cerana_QC*

# VEP
mkdir -p ${DUMP}/logs/VEP
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
  -o ${DUMP}/logs/VEP \
  -e ${DUMP}/logs/VEP \
  ${PIPE}/vep.sh \
    -i ${OUT}/Apis-cerana.vcf  \
    -o ${DUMP}/VEP/Apis-cerana.vcf 





bcftools view -H --min-ac 1:alt1 -G -O v -o ${OUT}/Apis-cerana-INFO.vcf \
  ${OUT}/Apis-cerana.vcf

sed 's/;/\t/g' ${OUT}/Apis-cerana-INFO.vcf | cut -f1-5,8,9 - | sed 's/AC=//g' - | sed 's/AN=//g' - | awk -F $'\t' 'BEGIN {OFS = FS} {print $1, $2, $3, $4, $5, $6/$7}' - > ${OUT}/Apis-cerana_AF.vcf

GTF=/save/dwragg/Apis/Apis_mellifera.GCA_000002195.1.25.gtf
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 1 \
  -o ${DUMP}/logs -e ${DUMP}/logs  \
  /work/dwragg/Analysis/cds/cds.sh \
    -e ${DUMP}/VEP/Apis-cerana.vcf_vep_output.txt  \
    -v ${OUT}/Apis-cerana_AF.vcf \
    -g ${GTF}.clean \
    -o ${DUMP}/cds/ACER



# Run locally
time ${DUMP}/cds/cds.py \
    -e ${DUMP}/VEP/Apis-cerana.vcf_vep_output.txt  \
    -v ${OUT}/Apis-cerana_AF.vcf \
    -g ${GTF}.clean \
    -o ${DUMP}/cds/ACER







chr=16
pos1=387001
pos2=793000
awk -v VAR1=${chr} -v VAR2=${pos1} -v VAR3=${pos2} '$1 ==VAR1 && $2 >=VAR2 && $2 <=VAR3' Apis-cerana.vcf





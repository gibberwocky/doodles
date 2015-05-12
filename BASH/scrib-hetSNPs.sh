#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a


# Load modules
module load bioinfo/bcftools

# Paths to main folders
DUMP=/home/dwragg/work/Analysis
PIPE=/save/seqapipop/Scripts
PERL=${PIPE}/Perl
SEQA=/save/seqapipop/Data/Apis-mellifera
AMEL=/work/dwragg/Apis-mellifera
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"

POP="JFM"

#mkdir -p ${DUMP}/hetSNPs/${POP}
#SAMPLES=($(cat ${DUMP}/vcfs/lists/${POP}_samples.list))

#for SAMPLE in ${SAMPLES[@]}
#do
#  mkdir -p ${DUMP}/hetSNPs/${POP}/${SAMPLE}
#  ln -fs ${SEQA}/SeqApiPop/${POP}/${SAMPLE}/${SAMPLE}_bootstrap* \
#    ${DUMP}/hetSNPs/${POP}/${SAMPLE}
#  mkdir -p ${DUMP}/logs/${SAMPLE}
#  cd ${DUMP}
#  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
#    -o ${DUMP}/logs/${SAMPLE} \
#    -e ${DUMP}/logs/${SAMPLE} \
#    ${PIPE}/snps-diploid.sh -s ${SAMPLE} -f ${DUMP}/hetSNPs/${POP} \
#      -o ${DUMP}/hetSNPs/${POP}
#done

# Merge samples and retain sites with one or more het genotypes
VCFS=(${DUMP}/hetSNPs/*/*/vcfs/*_clean.vcf.gz)
bcftools merge -m both ${VCFS[@]} \
  | bcftools view -m2 -M2 --types snps \
  --include 'DP>=9' \
  -g het -O z -o ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz
tabix -fp vcf ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz

# Reduce to chromosomes 1-16
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3' \
  -o ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf \
  ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz
bgzip -f ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf
tabix -fp vcf ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz

# Remove unnecessary annotations
bcftools annotate \
  -x FILTER,INF/Dels,INF/FS,INF/HaplotypeScore,INF/MQ0F,INF/MQSB,INF/MQ,\
INF/QD,INF/SF,INF/SOR,INF/BaseQRankSum,INF/MQRankSum,INF/SGB,INF/VDB,\
INF/ReadPosRankSum,INF/AF,INF/MLEAC,INF/MLEAF,INF/MQ0,INF/DP4,\
INF/BQB,INF/HOB,INF/ICB,INF/MQB,INF/RPB,FMT/GQ,FMT/FT,FMT/PL,FMT/AD,\
FMT/GL,FMT/GOF,FMT/NR,FMT/NV \
  -O v -o ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf \
  ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz
bgzip -f ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf
tabix -fp vcf ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz

# Remove duplicate SNP entries
bcftools norm -D -O v -o ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf \
  ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz 
bgzip -f ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf
tabix -fp vcf ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz

# Merge multialleleics
bcftools norm -m +any -O v -o ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf \
  ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz 
rm ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz 
rm ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf.gz.tbi

# Change from accession to chromosome naming and fix _ in sample names
${PIPE}/accession-to-chr.sh -i ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf
sed -i '/^#CHROM/s/_/-/g' ${DUMP}/hetSNPs/JFMOTH-hetSNPs.vcf


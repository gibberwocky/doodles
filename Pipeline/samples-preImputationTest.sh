#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

# Paths to main folders
DUMP=/home/dwragg/work/Analysis
PIPE=/save/seqapipop/Scripts
PERL=${PIPE}/Perl
SEQA=/save/seqapipop/Data/Apis-mellifera
AMEL=/work/dwragg/Apis-mellifera

# Merged master sites: 		HARP-JFM-OTH.vcf.gz
# Harpur_preImp.vcf.gz wc -l = 10950024 (header = 5338)
# JFMOTH_preImp.vcf.gz wc -l = 10981674 (header = 5338)

# This is combined Harpur and JFM-OTH SNPs post-filtered for depth of coverage

# Output path and master site file
OUT=${DUMP}/vcfs/HARP-JFM-OTH-15012015
sites=${DUMP}/sites/HARP-JFM-OTH.vcf.gz

# As not imputing, can pool samples into 2 runs: (1) haploid and (2) diploid
Harpur=(${VCF_A[@]} ${VCF_M[@]} ${VCF_C[@]} ${VCF_O[@]})
printf "%s\n" "${Harpur[@]}" > ${DUMP}/vcfs/lists/Harpur_vcfs.list
Harpur=(${BAM_A[@]} ${BAM_M[@]} ${BAM_C[@]} ${BAM_O[@]})
printf "%s\n" "${Harpur[@]}" > ${DUMP}/vcfs/lists/Harpur_bams.list
cat ${DUMP}/vcfs/lists/JFM_vcfs.list ${DUMP}/vcfs/lists/OTH_vcfs.list \
  > ${DUMP}/vcfs/lists/JFMOTH_vcfs.list
cat ${DUMP}/vcfs/lists/JFM_bams.list ${DUMP}/vcfs/lists/OTH_bams.list \
  > ${DUMP}/vcfs/lists/JFMOTH_bams.list

# Run (pay attention to ploidy flag!)
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/geno-sites.sh -o ${OUT} -s ${sites} -p 2 \
      -b ${DUMP}/vcfs/lists/Harpur_bams.list \
      -n Harpur


# ==============================================================================
# Create diploid drones using JFMOTH_preImp_chrs16_qc file
# ==============================================================================

OUT=${DUMP}/vcfs/HARP-JFM-OTH-15012015
RAW_VCF="JFMOTH_preImp"

# Reduce to chromosomes 1-16
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3' \
  -o ${OUT}/${RAW_VCF}_chrs16.vcf \
  ${OUT}/${RAW_VCF}.vcf.gz

# Change from accession to chromosome naming and fix _ in sample names
${PIPE}/accession-to-chr.sh -i ${OUT}/${RAW_VCF}_chrs16.vcf 
sed -i '/^#CHROM/s/_/-/g' ${OUT}/${RAW_VCF}_chrs16.vcf 

# Perform QC in Plink
# Slightly different to the main QC. purpose here is to remove low
# calling genotypes, not to worry about MAF because an allele can be fixed
# in a population
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
${PLINK}/plink --vcf ${OUT}/${RAW_VCF}_chrs16.vcf \
  --keep-allele-order \
  --a2-allele ${OUT}/${RAW_VCF}_chrs16.vcf 4 3 '#' \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --geno 0.1 \
  --mind 1 \
  --out ${OUT}/${RAW_VCF}_chrs16_qc \
  --recode vcf-iid

bgzip -f ${OUT}/${RAW_VCF}_chrs16_qc.vcf
tabix -fp vcf ${OUT}/${RAW_VCF}_chrs16_qc.vcf.gz

# Retain log
mv ${OUT}/${RAW_VCF}_chrs16_qc.log ${OUT}/${RAW_VCF}_chrs16_qc_for_diploids.log

# Remove unnecessary files
rm ${OUT}/${RAW_VCF}_chrs16_qc.nosex

# Impute missing 5% genotypes
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 4 \
  -o ${DUMP}/logs \
  -e ${DUMP}/logs \
  ${PIPE}/impute.sh \
    -v ${OUT}/${RAW_VCF}_chrs16_qc.vcf.gz \
    -o ${OUT}/${RAW_VCF}_chrs16_qc_beagle

# Once imputed, add contig lines to VCF header if not already present
gunzip ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf.gz 
sed -i '/^#C.*$/i \!contigs' ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
sed -i '/^!contigs.*/r /save/seqapipop/Analysis/sites/VCF_contigs.header' \
  ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
sed -i '/^!contigs.*/d' ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
bgzip -f ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
tabix -fp vcf ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf.gz

# Identify samples for diploidy
POP="OTH"
SUFFIX="_preImp"
readarray HAPS < ${DUMP}/vcfs/lists/${POP}_samples_fixed.list

# Create file to record sample IDs for merging
rm -f ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
touch ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
# Loop through all pairs
x=0
for ((i=0; i<${#HAPS[@]}; i=i+2))
do
  # Write sample ID to file
  printf "${OUT}/%s.vcf.gz\n" "${POP}-DIP${x}" \
    >> ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
  # Create directory for logs
  mkdir -p ${DUMP}/logs/${POP}-DIP${x}
  # Create diploid from samples i and i+1
  qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/${POP}-DIP${x} \
    -e ${DUMP}/logs/${POP}-DIP${x} \
    ${PIPE}/diploid-drone.sh \
      -v ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf.gz \
      -a ${HAPS[${i}]} \
      -b ${HAPS[${i}+1]} \
      -s ${POP}-DIP${x} \
      -o ${OUT}
  x=$((x+1))
done

# Remove dodgy diploids ( 0|-  1|-  -|0  -|1 ) = those with AN!=2
cd ${OUT}
POP="JFM"
ID=
ID=(${POP}-DIP*.vcf.gz)
ID=(${ID[@]/%.vcf.gz})
for VCF in ${ID[@]}
do
  bcftools view --type snps --exclude 'AN!=2' -O z -o ${VCF}-clean.vcf.gz \
    ${VCF}.vcf.gz
  bcftools annotate -x FORMAT -O z -o ${VCF}.vcf.gz \
    ${VCF}-clean.vcf.gz
  tabix -fp vcf ${VCF}.vcf.gz
  rm ${VCF}-clean.vcf.gz
done

# Merge diploids into single VCF
bcftools merge -m snps -O z \
  -l ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list \
  -o ${OUT}/${POP}_diploid.vcf.gz
tabix -fp vcf ${OUT}/${POP}_diploid.vcf.gz

# Remove unnecessary files
readarray FILES < ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
rm ${FILES[@]}
tmp=$(printf "%s.tbi " ${FILES[@]})
rm ${tmp}


# ==============================================================================
# Merge datasets
# ==============================================================================
# Diploid dataset for Admixture, etc.
bcftools merge -m snps -O z \
  ${OUT}/Harpur_preImp.vcf.gz \
  ${OUT}/JFM_diploid.vcf.gz \
  ${OUT}/OTH_diploid.vcf.gz \
  -o ${OUT}/Harpur-dipJFMOTH_preImp.vcf.gz
tabix -fp vcf ${OUT}/Harpur-dipJFMOTH_preImp.vcf.gz

# Haploid dataset for Admixture to compare against diploid dataset
bcftools merge -m snps -O z \
  ${OUT}/Harpur_preImp.vcf.gz \
  ${OUT}/JFMOTH_preImp.vcf.gz \
  -o ${OUT}/Harpur-hapJFMOTH_preImp.vcf.gz
tabix -fp vcf ${OUT}/Harpur-hapJFMOTH_preImp.vcf.gz


# ==============================================================================
# QC Filter 
# ==============================================================================

OUT=${DUMP}/vcfs/HARP-JFM-OTH-15012015
RAW_VCF="JFMOTH_preImp"

# Reduce to chromosomes 1-16
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3' \
  -o ${OUT}/${RAW_VCF}_chrs16.vcf \
  ${OUT}/${RAW_VCF}.vcf.gz

# Change from accession to chromosome naming and fix _ in sample names
${PIPE}/accession-to-chr.sh -i ${OUT}/${RAW_VCF}_chrs16.vcf 
sed -i '/^#CHROM/s/_/-/g' ${OUT}/${RAW_VCF}_chrs16.vcf 

# Perform QC in Plink
# CAUTION, --geno and --mind are exclusion threholds, eg:
# --geno 0.1 will remove all SNPs with less than 10% call rate
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
${PLINK}/plink --vcf ${OUT}/${RAW_VCF}_chrs16.vcf \
  --keep-allele-order \
  --a2-allele ${OUT}/${RAW_VCF}_chrs16.vcf 4 3 '#' \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --maf 0.05 \
  --geno 0.1 \
  --mind 1 \
  --out ${OUT}/${RAW_VCF}_chrs16_qc \
  --recode vcf-iid

bgzip -f ${OUT}/${RAW_VCF}_chrs16_qc.vcf
tabix -fp vcf ${OUT}/${RAW_VCF}_chrs16_qc.vcf.gz

# Retain log
mv ${OUT}/${RAW_VCF}_chrs16_qc.log ${OUT}/${RAW_VCF}_chrs16_qc_initial.log

# Remove unnecessary files
rm ${OUT}/${RAW_VCF}_chrs16_qc.nosex

# Impute missing 5% genotypes
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 4 \
  -o ${DUMP}/logs \
  -e ${DUMP}/logs \
  ${PIPE}/impute.sh \
    -v ${OUT}/${RAW_VCF}_chrs16_qc.vcf.gz \
    -o ${OUT}/${RAW_VCF}_chrs16_qc_beagle

# Once imputed, add contig lines to VCF header if not already present
gunzip ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf.gz 
sed -i '/^#C.*$/i \!contigs' ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
sed -i '/^!contigs.*/r /save/seqapipop/Analysis/sites/VCF_contigs.header' \
  ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
sed -i '/^!contigs.*/d' ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
bgzip -f ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
tabix -fp vcf ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf.gz



# ==============================================================================
# Run analyses (analysis-preImputationTests.sh)
# ==============================================================================







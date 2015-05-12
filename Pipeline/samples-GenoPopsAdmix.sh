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
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"


# Update lists
# Refer to paths.sh to update BAM and VCF lists

# Output path and master site file
sites=${DUMP}/sites/HARP-JFM-OTH.vcf.gz
OUT=${DUMP}/vcfs/HARP-JFM-OTH-BR-COR-OUE-BER-SAR-SAV
mkdir -p ${OUT}

# Combine target population bam lists into a single file
cat ${DUMP}/vcfs/lists/BR_bams.list ${DUMP}/vcfs/lists/BER_bams.list \
  ${DUMP}/vcfs/lists/COR_bams.list ${DUMP}/vcfs/lists/OUE_bams.list \
  ${DUMP}/vcfs/lists/SAR_bams.list ${DUMP}/vcfs/lists/SAV_bams.list \
  > ${DUMP}/vcfs/lists/BR-BER-COR-OUE-SAR-SAV_bams.list

# Run (pay attention to ploidy flag!)
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/geno-sites.sh -o ${OUT} -s ${sites} -p 1 \
      -b ${DUMP}/vcfs/lists/BR-BER-COR-OUE-SAR-SAV_bams.list \
      -n SAR-SAV

# Haploid dataset for Admixture to compare against diploid dataset
bcftools merge -m snps -O z \
  ${OUT}/BR-BER-COR-OUE_preImp.vcf.gz \
  ${OUT}/SAR-SAV_preImp.vcf.gz \
  -o ${OUT}/BR-BER-COR-OUE-SAR-SAV_preImp.vcf.gz
tabix -fp vcf ${OUT}/BR-BER-COR-OUE-SAR-SAV_preImp.vcf.gz



# ===================================================================
# Create diploid drones
# ===================================================================

RAW_VCF="BR-BER-COR-OUE-SAR-SAV_preImp"

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
cp ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf \
  ${OUT}/${RAW_VCF}_chrs16_qc_beagle_backup.vcf
sed -i '/^#C.*$/i \!contigs' ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
sed -i '/^!contigs.*/r /save/seqapipop/Analysis/sites/VCF_contigs.header' \
  ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
sed -i '/^!contigs.*/d' ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
bgzip -f ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf 
tabix -fp vcf ${OUT}/${RAW_VCF}_chrs16_qc_beagle.vcf.gz



# Identify samples for diploidy
# COR[37], OUE[40], BER[19], SAR[18], SAV[31], BR[2]
POP="COR"

# Generate sample lists (necessary due to different number of samples now
# available for some of the populations, since imputation)
#TMP=($(cat ${DUMP}/vcfs/lists/${POP}_bams.list))
#TMP=(${TMP[@]%_bootstrap.bam})
#TMP=(${TMP[@]#${AMEL}/SeqApiPop/BR/bams/})
#printf "%s\n" "${TMP[@]}" > ${DUMP}/vcfs/lists/${POP}_samples.list
#sed 's/_/-/g' ${DUMP}/vcfs/lists/${POP}_samples.list \
#  > ${DUMP}/vcfs/lists/${POP}_samples_fixed.list

SUFFIX="_preImp"
readarray HAPS < ${DUMP}/vcfs/lists/${POP}_samples_fixed.list

# Make sure HAPS array is an even size
HAPL=$(( ${#HAPS[@]} % 2 ))
if [ $HAPL -eq 0 ]; then
  echo " even ${#HAPS[@]}"
else
  echo " odd ${#HAPS[@]}"
  HAPS=(${HAPS[@]:0:$((${#HAPS[@]}-1))}) 
fi

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
POP="OUE"
ID=
ID=(${POP}-DIP*.vcf.gz)
ID=(${ID[@]/%.vcf.gz})
for VCF in ${ID[@]}
do
  # Create diploid from samples i and i+1
  qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/${VCF} \
    -e ${DUMP}/logs/${VCF} \
    ${PIPE}/dodgydips.sh -o ${OUT} -v ${VCF}
done

# Merge diploids into single VCF and remove unnecessary files
POP="OUE"
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
  -o ${DUMP}/logs -e ${DUMP}/logs \
  ${PIPE}/mergeVCFs.sh -o ${OUT} -p ${POP} -i ${DUMP}/vcfs/lists

# =================   end of diploid construction   =================


# Diploid dataset for Admixture, etc.
OUTa=~/work/Analysis/vcfs/HARP-JFM-OTH-15012015
bcftools merge -m snps -O z \
  ${OUTa}/Harpur_preImp_chrs16_qc.vcf.gz \
  ${OUTa}/JFM_diploid.vcf.gz \
  ${OUTa}/OTH_diploid.vcf.gz \
  ${OUT}/BR_diploid.vcf.gz \
  ${OUT}/BER_diploid.vcf.gz \
  ${OUT}/COR_diploid.vcf.gz \
  ${OUT}/OUE_diploid.vcf.gz \
  ${OUT}/SAR_diploid.vcf.gz \
  ${OUT}/SAV_diploid.vcf.gz \
  -o ${OUT}/Mixedbag.vcf.gz
tabix -fp vcf ${OUT}/Mixedbag.vcf.gz



# Admixture analysis

VCF="Mixedbag"
gunzip ${OUT}/${VCF}.vcf.gz
mkdir -p ${DUMP}/vcfs/plink/${VCF}

# Make BED file
${PLINK}/plink --vcf ${OUT}/${VCF}.vcf \
  --keep-allele-order \
  --a2-allele ${OUT}/${VCF}.vcf 4 3 '#' \
  --allow-no-sex \
  --allow-extra-chr \
  --mind 1 \
  --geno 0.1 \
  --maf 0.05 \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/plink/${VCF}/${VCF} \
  --make-bed

# --a2-allele: 10319963 assignments made.
# 1402502 variants removed due to missing genotype data (--geno).
# 6303979 variants removed due to MAF threshold(s) (--maf/--max-maf).
# 2613482 variants and 134 samples pass filters and QC.


# Run Admixture on cluster
k=7 # Min K
K=8 # Max K
N=1 # Number of bootstraps
SUFFIX=""
mkdir -p ${DUMP}/vcfs/admixture/${VCF}${SUFFIX}
for ((i=k; i<$((${K}+1)); i=i+1))
do
  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/admixture.sh -i ${DUMP}/vcfs/plink/${VCF}/${VCF}${SUFFIX} \
      -k ${i} -K ${i} -n ${N} -o ${DUMP}/vcfs/admixture/${VCF}${SUFFIX}/ \
      -s ${VCF}${SUFFIX}
done


# Copy fam file to same folder incase needed for plotting
cp ${DUMP}/vcfs/plink/${VCF}/${VCF}${SUFFIX}.fam \
  ${DUMP}/vcfs/admixture/${VCF}${SUFFIX}

# Identify optimal K
grep -h CV ${DUMP}/vcfs/admixture/${VCF}${SUFFIX}/${VCF}-*.out




# =============================================================================
# Plink to R - for phylogenetic analysis
# =============================================================================

# This assumes that the Plink file used in the ADMIXTURE analysis is intacted
# and in the ${DUMP}/vcfs/plink folder, and that the ADMIXTURE results have
# been transferred to their own folder in ${DUMP}/vcfs/admixture, including
# the plink fam file.

SITES=${OUT}
VCF="Mixedbag"
INDIR=${DUMP}/vcfs/plink/${VCF}

# FAM to PHENO (id	sex	 phene)
cut -d' ' -f1,3,4 ${INDIR}/${VCF}.fam > ${INDIR}/${VCF}-R.tmp
header=('id' 'sex' 'pop')
echo ${header[@]} | cat - ${INDIR}/${VCF}-R.tmp > ${INDIR}/${VCF}-R.pheno
rm ${INDIR}/${VCF}-R.tmp

# Format Plink for GenABEL import
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${VCF}/${VCF} \
  --keep-allele-order \
  --allow-no-sex \
  --allow-extra-chr \
  --mind 1 \
  --geno 0.1 \
  --maf 0.05 \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --recode 12 \
  --out ${DUMP}/vcfs/plink/${VCF}/${VCF}-R

# Remove genetic distance from map file
cut -f1,2,4 ${DUMP}/vcfs/plink/${VCF}/${VCF}-R.map \
  > ${DUMP}/vcfs/plink/${VCF}/${VCF}-R.tmp
mv ${DUMP}/vcfs/plink/${VCF}/${VCF}-R.tmp ${DUMP}/vcfs/plink/${VCF}/${VCF}-R.map

# Proceed in R




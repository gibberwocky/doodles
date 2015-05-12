#!/bin/bash
# IEKNGHER1ujm

# BER [17*]	BR [7]		Corse [37]	JFM [30]
# OTH [30]	Ouessant [39]

# Paths to main folders
DUMP="/home/dwragg/work/Analysis"
PIPE=/save/seqapipop/Scripts
PERL=${PIPE}/Perl

# Path to linked data
SEQA=/save/seqapipop/Data/Apis-mellifera
#SEQA=/genphyse/cytogen/seqapipop/Data
AMEL=/work/dwragg/Apis-mellifera

# Soft link data
NAME="SeqApiPop/Sarthe"
mkdir -p ${AMEL}/${NAME}/vcfs
mkdir -p ${AMEL}/${NAME}/bams
ln -fs ${SEQA}/${NAME}/*/vcfs/*_clean.vcf.gz ${AMEL}/${NAME}/vcfs/
ln -fs ${SEQA}/${NAME}/*/*_bootstrap* ${AMEL}/${NAME}/bams/
for ID in ${AMEL}/${NAME}/vcfs/*.vcf.gz; do 
  echo "${ID}"
  tabix -fp vcf ${ID}
done
# Sometimes the tabix command needs running independently to avoid error
# linked to consecutive blocks being out of sequence:
cd ${AMEL}/${NAME}/vcfs
ID=(*)
for SAMPLE in ${ID[@]}
do
  tabix -fp vcf ${SAMPLE}
done

# Links to sites
ln -fs /save/seqapipop/Analysis/sites/* ${DUMP}/sites

# Paths to linked VCF data from sequence read archive
VCF_HARP=(${AMEL}/Harpur/vcfs/*_clean.vcf.gz)
VCF_WALL=(${AMEL}/Wallberg/vcfs/*_clean.vcf.gz)
VCF_LIU=(${AMEL}/Liu/vcfs/*_clean.vcf.gz)

# Paths to linked BAMs from sequence read archive
BAM_HARP=(${AMEL}/Harpur/bams/*_bootstrap.bam)
BAM_WALL=(${AMEL}/Wallberg/bams/*_bootstrap.bam)
BAM_LIU=(${AMEL}/Liu/bams/*_bootstrap.bam)

# Paths to linked VCF data from SeqApiPop
VCF_JFM=(${AMEL}/SeqApiPop/JFM/vcfs/*_clean.vcf.gz)
VCF_OTH=(${AMEL}/SeqApiPop/OTH/vcfs/*_clean.vcf.gz)
VCF_BR=(${AMEL}/SeqApiPop/BR/vcfs/*_clean.vcf.gz)
VCF_COR=(${AMEL}/SeqApiPop/Corse/vcfs/*_clean.vcf.gz)
VCF_OUE=(${AMEL}/SeqApiPop/Ouessant/vcfs/*_clean.vcf.gz)
VCF_BER=(${AMEL}/SeqApiPop/Berlin/vcfs/*_clean.vcf.gz)
VCF_SAR=(${AMEL}/SeqApiPop/Sarthe/vcfs/*_clean.vcf.gz)
VCF_SAV=(${AMEL}/SeqApiPop/Savoie/vcfs/*_clean.vcf.gz)

# Paths to linked BAMs from SeqApiPop
BAM_JFM=(${AMEL}/SeqApiPop/JFM/bams/*_bootstrap.bam)
BAM_OTH=(${AMEL}/SeqApiPop/OTH/bams/*_bootstrap.bam)
BAM_BR=(${AMEL}/SeqApiPop/BR/bams/*_bootstrap.bam)
BAM_COR=(${AMEL}/SeqApiPop/Corse/bams/*_bootstrap.bam)
BAM_OUE=(${AMEL}/SeqApiPop/Ouessant/bams/*_bootstrap.bam)
BAM_BER=(${AMEL}/SeqApiPop/Berlin/bams/*_bootstrap.bam)
BAM_SAR=(${AMEL}/SeqApiPop/Sarthe/bams/*_bootstrap.bam)
BAM_SAV=(${AMEL}/SeqApiPop/Savoie/bams/*_bootstrap.bam)

# Write VCF and BAM lists
printf "%s\n" "${VCF_HARP[@]}" > ${DUMP}/vcfs/lists/Harpur_vcfs.list
printf "%s\n" "${BAM_HARP[@]}" > ${DUMP}/vcfs/lists/Harpur_bams.list
printf "%s\n" "${VCF_WALL[@]}" > ${DUMP}/vcfs/lists/Wallberg_vcfs.list
printf "%s\n" "${BAM_WALL[@]}" > ${DUMP}/vcfs/lists/Wallberg_bams.list
printf "%s\n" "${VCF_LIU[@]}" > ${DUMP}/vcfs/lists/Liu_vcfs.list
printf "%s\n" "${BAM_LIU[@]}" > ${DUMP}/vcfs/lists/Liu_bams.list
printf "%s\n" "${VCF_JFM[@]}" > ${DUMP}/vcfs/lists/JFM_vcfs.list
printf "%s\n" "${BAM_JFM[@]}" > ${DUMP}/vcfs/lists/JFM_bams.list
printf "%s\n" "${VCF_OTH[@]}" > ${DUMP}/vcfs/lists/OTH_vcfs.list
printf "%s\n" "${BAM_OTH[@]}" > ${DUMP}/vcfs/lists/OTH_bams.list
printf "%s\n" "${VCF_BR[@]}" > ${DUMP}/vcfs/lists/BR_vcfs.list
printf "%s\n" "${BAM_BR[@]}" > ${DUMP}/vcfs/lists/BR_bams.list
printf "%s\n" "${VCF_COR[@]}" > ${DUMP}/vcfs/lists/COR_vcfs.list
printf "%s\n" "${BAM_COR[@]}" > ${DUMP}/vcfs/lists/COR_bams.list
printf "%s\n" "${VCF_OUE[@]}" > ${DUMP}/vcfs/lists/OUE_vcfs.list
printf "%s\n" "${BAM_OUE[@]}" > ${DUMP}/vcfs/lists/OUE_bams.list
printf "%s\n" "${VCF_BER[@]}" > ${DUMP}/vcfs/lists/BER_vcfs.list
printf "%s\n" "${BAM_BER[@]}" > ${DUMP}/vcfs/lists/BER_bams.list
printf "%s\n" "${VCF_SAR[@]}" > ${DUMP}/vcfs/lists/SAR_vcfs.list
printf "%s\n" "${BAM_SAR[@]}" > ${DUMP}/vcfs/lists/SAR_bams.list
printf "%s\n" "${VCF_SAV[@]}" > ${DUMP}/vcfs/lists/SAV_vcfs.list
printf "%s\n" "${BAM_SAV[@]}" > ${DUMP}/vcfs/lists/SAV_bams.list

# Write samples list
NAME="SeqApiPop/Ouessant"
OUTNOM="OUE"
cd ${SEQA}/${NAME}
ID=(*)
printf "%s\n" "${ID[@]}" > ${DUMP}/vcfs/lists/${OUTNOM}_samples.list
sed 's/_/-/g' ${DUMP}/vcfs/lists/${OUTNOM}_samples.list \
  > ${DUMP}/vcfs/lists/${OUTNOM}_samples_fixed.list





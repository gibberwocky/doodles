#!/bin/bash
# IEKNGHER1ujm


# ==============================================================================
# Check results carefully!
# Flagstat output changed between mapping JFM/OTH/BR and COR/OUE
# Perl script was recording the wrong value for single/paired read counts
# ==============================================================================


# Set-up
DUMP="/home/dwragg/work/Analysis"
PIPE="/home/dwragg/work/Pipeline"
SCRIPTS="/save/seqapipop/Scripts"

# Pre-reqs
POP="Corse"
INDIR="/save/seqapipop/Data/Apis-mellifera/SeqApiPop/${POP}"
#INDIR="/save/seqapipop/Data/Apis-mellifera/${POP}"

cd ${INDIR}
SAMPLES=(*)

cd ${DUMP}
printf "%s\n" "${SAMPLES[@]}" > ${DUMP}/tmp


# SeqApiPop data
perl ${SCRIPTS}/Perl/StatSummary.pl ${INDIR} ${DUMP}/tmp > ${DUMP}/${POP}.stats

# No insert mertrics available for the Harpur data (not paired read data)
#perl ${SCRIPTS}/Perl/StatSummary-Harpur.pl ${INDIR} ${DUMP}/tmp > ${DUMP}/${POP}.stats

# Remove temporary file
rm tmp




# Calculate total SNPs per group (SeqApiPop vs Harpur)
AMEL=/work/dwragg/Apis-mellifera
VCF_HARP=(${AMEL}/Harpur/vcfs/*_clean.vcf.gz)
VCF_JFM=(${AMEL}/SeqApiPop/JFM/vcfs/*_clean.vcf.gz)
VCF_OTH=(${AMEL}/SeqApiPop/OTH/vcfs/*_clean.vcf.gz)


# N*3*average coverage
SAMPLES=(${VCF_JFM[@]} ${VCF_OTH[@]})
bcftools merge ${SAMPLES[@]} \
  | bcftools view --types snps \
  --min-af 0.1 --include 'DP>=60 & DP<=1278' \
  -G -M 2 -O z -o ${DUMP}/JFM-OTH-n3d.vcf.gz -
tabix -p vcf ${DUMP}/JFM-OTH-n3d.vcf.gz

# Copy to /save/seqapipop, remove source and link
cp ${DUMP}/Harpur-n3d.vcf* /save/seqapipop/Analysis/sites



ln -fs /save/seqapipop/Analysis/sites/* ${DUMP}/sites


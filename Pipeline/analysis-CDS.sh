#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

PIPE=/save/seqapipop/Scripts
DUMP=/work/dwragg/Analysis
SITES=${DUMP}/vcfs/HARP-JFM-OTH-15012015
VCF="JFMOTH_preImp_chrs16_qc_beagle"
VCFT="/usr/local/bioinfo/src/vcftools/vcftools_0.1.12a/bin"


# =============================================================================
# Generate INFO VCF files for each subset of samples per population, where
# at least 1 copy of the alternative allele is present. This is to ensure
# that the SNP was called at least once in the population. The sed part of
# the script re-calculates the AF field based on the data subset AC/AN fields.
# This is necessary becasuse even after sub-setting the data, the original
# AF field contuinues to relate to the data prior sub-setting.
# =============================================================================
POPS=('JFM' 'OTH')
for POP in ${POPS[@]}
do
 # Output population-specific INFO (AC and AN only; AF relates to whole data)
  bcftools view -H --min-ac 1:alt1 -G -O v -o ${SITES}/${VCF}_${POP}_info.vcf \
    -S ${DUMP}/vcfs/lists/${POP}_samples_fixed.list \
    ${SITES}/${VCF}.vcf.gz
  # Reduce down to just chromosome annotation and re-calculated AF (AC/AN)
  sed 's/;/\t/g' ${SITES}/${VCF}_${POP}_info.vcf | cut -f1-5,11-12 - | sed 's/AC=//g' - | sed 's/AN=//g' - | awk -F $'\t' 'BEGIN {OFS = FS} {print $1,  $2, $3, $4, $5, $6/$7}' - > ${SITES}/${VCF}_${POP}_AF.vcf
  # Remove surplous file
  rm ${SITES}/${VCF}_${POP}_info.vcf
done

# =============================================================================
# Run VEP on the results for each inidividual population
# =============================================================================
POPS=('JFM' 'OTH')
for POP in ${POPS[@]}
do
  mkdir -p ${DUMP}/logs/VEP
  qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/VEP \
    -e ${DUMP}/logs/VEP \
    ${PIPE}/vep.sh \
      -i ${SITES}/${VCF}_${POP}_AF.vcf \
      -o ${DUMP}/VEP/${VCF}_${POP}
done

# =============================================================================
# Parse GTF file to strip out unnecessary columns and gene id
# =============================================================================
GTF=/save/dwragg/Apis/Apis_mellifera.GCA_000002195.1.25.gtf
awk -F'\t' 'BEGIN {OFS = FS}; {split($9, g, ";"); printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $7, g[1]);}' ${GTF} | sed 's/gene_id//g' - | sed 's/\"//g' - | sed 's/\t */\t/g' - > ${GTF}.clean

# =============================================================================
# Run CDS anlaysis to compute alpha per gene per population
# =============================================================================
GTF=/save/dwragg/Apis/Apis_mellifera.GCA_000002195.1.25.gtf
POPS=('JFM' 'OTH')
for POP in ${POPS[@]}
do
  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 1 \
    -o ${DUMP}/logs -e ${DUMP}/logs  \
    ${PIPE}/cds.sh \
      -e ${DUMP}/VEP/${VCF}_${POP}_vep_output.txt  \
      -v ${SITES}/${VCF}_${POP}_AF.vcf \
      -g ${GTF}.clean \
      -o ${DUMP}/cds/${POP}
done

# =============================================================================
# Recover non-synonmous (N) and synonmous (S) SNPs in haplotype format
# For plotting haplotypes at genes of interest
# =============================================================================
POPS=('JFM' 'OTH')
for POP in ${POPS[@]}
do
  awk -F'\t' 'BEGIN {OFS = FS}; { print $14; }' ${DUMP}/cds/${POP}.csv \
    | sed 's/,/\n/g' - | sed '/^\s*$/d' - > ${DUMP}/cds/${POP}-N.snps
  awk -F'\t' 'BEGIN {OFS = FS}; { print $15; }' ${DUMP}/cds/${POP}.csv \
    | sed 's/,/\n/g' - | sed '/^\s*$/d' - > ${DUMP}/cds/${POP}-S.snps
  ${VCFT}/vcftools --gzvcf ${SITES}/${VCF}.vcf.gz \
    --snps ${DUMP}/cds/${POP}-N.snps \
    --012 --out ${DUMP}/cds/${POP}-N.snps
  ${VCFT}/vcftools --gzvcf ${SITES}/${VCF}.vcf.gz \
    --snps ${DUMP}/cds/${POP}-S.snps \
    --012 --out ${DUMP}/cds/${POP}-S.snps
done














size=1000

# To subset data if needed for testing
head --l ${size} ${DUMP}/VEP/JFMOTH_preImp_chrs16_qc_beagle.vcf_vep_output.txt > ${DUMP}/cds/VEPh.txt
tail --l ${size} ${DUMP}/VEP/JFMOTH_preImp_chrs16_qc_beagle.vcf_vep_output.txt > ${DUMP}/cds/VEPt.txt
cat ${DUMP}/cds/VEPh.txt ${DUMP}/cds/VEPt.txt > ${DUMP}/cds/VEP.txt
rm ${DUMP}/cds/VEPh.txt
rm ${DUMP}/cds/VEPt.txt

head --l ${size} ${SITES}/JFMOTH_preImp_chrs16_qc_beagle_JFM_AF.vcf > ${DUMP}/cds/VCFh.txt
tail --l ${size} ${SITES}/JFMOTH_preImp_chrs16_qc_beagle_JFM_AF.vcf > ${DUMP}/cds/VCFt.txt
cat ${DUMP}/cds/VCFh.txt ${DUMP}/cds/VCFt.txt > ${DUMP}/cds/VCF.txt
rm ${DUMP}/cds/VCFh.txt
rm ${DUMP}/cds/VCFt.txt



vcffile="/work/dwragg/Analysis/cds/VCF.txt"
vepfile="/work/dwragg/Analysis/cds/VEP.txt"
gtffile="/save/dwragg/Apis/Apis_mellifera.GCA_000002195.1.25.gtf.clean"
outfile="/work/dwragg/Analysis/cds/test"

# Run locally
time ${PIPE}/cds.py \
    -e ${vepfile} \
    -v ${vcffile} \
    -g /save/dwragg/Apis/Apis_mellifera.GCA_000002195.1.25.gtf.clean \
    -o ${DUMP}/cds/test



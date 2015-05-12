#!/bin/bash
# IEKNGHER1ujm


# Set-up
module load bioinfo/bcftools
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PERL=${PIPE}/Perl
VCFT="/usr/local/bioinfo/src/vcftools/vcftools_0.1.12a/bin"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"

# Pre-imputation haploid JFM OTH data
SITES=${DUMP}/vcfs/HARP-JFM-OTH-15012015
VCF="Harpur-hapJFMOTH_preImp_chrs16_qc"

# Extract target population
POPS=('A_Harpur' 'M_Harpur' 'C_Harpur' 'O_Harpur' 'JFM' 'OTH')

for POP in ${POPS[@]}
do
  
  # Check which population to reference correct location for samples list
  if [[ "${POP}" == "JFM" || "${POP}" == "OTH" ]]
  then
    SAMPLES=${DUMP}/selection/hp/${POP}_samples.list
  else
    SAMPLES=${DUMP}/vcfs/lists/${POP}_samples.list
  fi

  # Extract Harpur data
  bcftools view --types snps -M 2 -O v \
    -o ${DUMP}/selection/vcft-nuc/${POP}.vcf \
    -S ${SAMPLES} \
    ${SITES}/${VCF}.vcf.gz

  # Change from accession to chromosome naming
  ${PIPE}/accession-to-chr.sh -i ${DUMP}/selection/vcft-nuc/${POP}.vcf
  sed -i '/^#CHROM/s/_/-/g' ${DUMP}/selection/vcft-nuc/${POP}.vcf

done

# =============================================================================
# VCFtools
#	- windowed Nucleotide divergence statistics
# =============================================================================

# Calculate nucleotdie divergence (pi) in overlapping windows
POPS=('A_Harpur' 'M_Harpur' 'C_Harpur' 'O_Harpur' 'JFM' 'OTH')
for POP in ${POPS[@]}
do
  ${VCFT}/vcftools --vcf ${DUMP}/selection/vcft-nuc/${POP}.vcf \
    --window-pi 5000 \
    --window-pi-step 1000 \
    --out ${DUMP}/selection/vcft-nuc/${POP}
done


# Output summary stats
POPS=('A_Harpur' 'M_Harpur' 'C_Harpur' 'O_Harpur' 'JFM' 'OTH')
for POP in ${POPS[@]}
do
  awk 'NR == 1 { sum=0 } { sum+=$5;} END {printf "Average: %f\n", sum/NR}' \
    ${POP}.*.pi
done

# Original manuscript		New draft, preImp
#A    Average: 0.007679		Average: 0.006282
#M    Average: 0.003360		Average: 0.003129
#C    Average: 0.002740		Average: 0.002557
#O    Average: 0.003679		Average: 0.003413
#JFM  Average: 0.004583		Average: 0.003904
#OTH  Average: 0.003968		Average: 0.003470





#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

# Paths to main folders
DUMP=/home/dwragg/work/Analysis
PIPE=/save/seqapipop/Scripts
PERL=${PIPE}/Perl

REF=/save/dwragg/Apis/Apis_mellifera.fa

CHRS=('NC_007070.3' 'NC_007071.3' 'NC_007072.3' 'NC_007073.3' \
'NC_007074.3' 'NC_007075.3' 'NC_007076.3' 'NC_007077.3' 'NC_007078.3' \
'NC_007079.3' 'NC_007080.3' 'NC_007081.3' 'NC_007082.3' 'NC_007083.3' \
'NC_007084.3' 'NC_007085.3')

POP="JFM"
OUT=${DUMP}/dcov
BAMS=${DUMP}/vcfs/lists/${POP}_bams.list

for CHR in ${CHRS[@]}
do
  # Calculate per-base depth of coverage for all samples in population
  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${OUT} \
    -e ${OUT} \
    ${PIPE}/basecov.sh -c ${CHR} -q 20 -r ${REF} \
      -b ${BAMS} \
      -o ${OUT}/${POP}-${CHR}.vcf.gz
done

# Alternative method, per sample, by piping into bedtools
#BAM=/work/dwragg/Apis-mellifera/SeqApiPop/OTH/bams/FL11_GTGAAA_L004_bootstrap.bam
#samtools view -q 20 ${BAM} | bedtools genomecov -ibam - -d > test

# Strip out genotypes and save as VCF
i=0
for CHR in ${CHRS[@]}
do
  bcftools view --no-header -G -O v -o ${OUT}/${POP}-${CHR}.vcf \
    ${OUT}/${POP}-${CHR}.vcf.gz
  i=$((i+1))
  sed -i "s/${CHR}/${i}/" ${OUT}/${POP}-${CHR}.vcf
  mv ${OUT}/${POP}-${CHR}.vcf ${OUT}/${POP}-${i}.vcf
done

# Bin DP into windows
CHRS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16')
for CHR in ${CHRS[@]}
do
  perl ${PERL}/VCF-winDP.pl ${OUT}/${POP}-${CHR}.vcf ${CHR} 5000 1000
done

# Concatenate files for each population
DP=(${POP}*.dp)
cat ${DP[@]} > ${POP}.dp
sed -i '1s/^/chr\tstart\tend\tsum_DP\tavg_DP\n/' ${POP}.dp




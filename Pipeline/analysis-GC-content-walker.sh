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
OUT=${DUMP}/GC
# 1) extract single sample into a VCF file
# 2) generate alternate reference fasta
# 3) Calculate GC in windows per sample
# 4) Generate mean for each window for each population


#===============================================================================
# Run GC Walker
#===============================================================================
POP="OTH"
SAMPLES=($(cat ${DUMP}/vcfs/lists/${POP}_samples_fixed.list))

for SAMPLE in ${SAMPLES[@]}
do

  # GC walker
  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/GCwalker.sh \
      -v ${DUMP}/vcfs/HARP-JFM-OTH-15012015/JFMOTH_preImp_chrs16_qc_beagle.vcf.gz \
      -o ${DUMP}/GC -s ${SAMPLE} -w 5000 -i 1000

done


#===============================================================================
# Average the results across samples in each population
#===============================================================================
POP="OTH"
SAMPLES=($(cat ${DUMP}/vcfs/lists/${POP}_samples_fixed.list))

awk '{
  avg[$1,$2,$3]+=$4; 
  i[FNR]=$1 SUBSEP $2 SUBSEP $3
}
END{
  for(k in i) sorted_i[j++]=k+0; 
  n=asort(i,sorted_i); 
  for(j=1; j<=n; j++) 
    print i[j]"\t" avg[i[j]]/(ARGC-2)
}' SUBSEP='\t' ${SAMPLES[@]/%/.gc} > ${POP}.gc



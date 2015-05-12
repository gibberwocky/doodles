#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

# Paths to main folders
DUMP="/home/dwragg/work/Analysis"
PIPE="/home/dwragg/work/Pipeline"
PERL="/save/dwragg/Scripts/Perl"


NAME="HARP-WALL"

mkdir -p ${DUMP}/logs/${NAME}
cd ${DUMP}
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
  -o ${DUMP}/logs/${NAME} \
  -e ${DUMP}/logs/${NAME} \
  ${PIPE}/snps-pool.sh -i ${PIPE}/params -b ${DUMP}/vcfs/lists/${NAME}.list -n ${NAME} -l ${DUMP}/logs -o ${DUMP}/test





# Plan is to see if calling SNPs from pooled BAMs gives better results
# in the Harpur-Wallberg Admixture output...

# May be possible to perform admixture on the immediate output, rather than
# to re-genotype samples...

# This takes FAR too long, will eat all CPU time


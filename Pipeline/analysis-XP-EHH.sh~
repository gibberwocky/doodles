#!/bin/bash
# IEKNGHER1ujm

# Set-up
module load bioinfo/bcftools
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
SITES=${DUMP}/vcfs/HARP-JFM-OTH-07112014
VCF=${DUMP}/selection/ehh/S2_hapJFMOTH_beagle.vcf.gz


# ============================================================================= 
# Prep files for SelScan iHS/XP-EHH analysis
# ============================================================================= 

# NCBI chromosome list
CHRS=('NC_007070.3' 'NC_007071.3' 'NC_007072.3' 'NC_007073.3' \
  'NC_007074.3' 'NC_007075.3' 'NC_007076.3' 'NC_007077.3' 'NC_007078.3' \
  'NC_007079.3' 'NC_007080.3' 'NC_007081.3' 'NC_007082.3' 'NC_007083.3' \
  'NC_007084.3' 'NC_007085.3' 'NC_001566.1')

# Temporary output vcf name
FILE=JFM-OTH.vcf

# Create logs directory
mkdir -p ${DUMP}/selection/ehh/logs

# Population to operate on (only complete one population at a time)
POP="JFM"

# Column in VCF file where consecutive POP samples start
# JFM = 10, OTH = 40
START=10

# Number of samples in population
N=30

# Run for each chromosome on cluster
for CHR in ${CHRS[@]}
do
  qsub -q unlimitq -l mem=4G -l h_vmem=32G -pe parallel_smp 1 \
    -o ${DUMP}/selection/ehh/logs \
    -e ${DUMP}/selection/ehh/logs \
    ${PIPE}/iHH-prep.sh \
      -p ${DUMP}/selection/ehh \
      -c ${CHR} \
      -o ${FILE} \
      -v ${VCF} \
      -i ${POP} \
      -x ${START} \
      -n ${N}
done



# ============================================================================= 
# SelScan
# https://github.com/szpiech/selscan
# ============================================================================= 
# Apis mellifera = ~ 28.7 kb per cM



# -----------------------------------------------------------------------------
# XP-EHH
# -----------------------------------------------------------------------------

POP1="JFM"
POP2="OTH"
mkdir -p ${DUMP}/logs/selscan
for CHR in ${CHRS[@]}
do
  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs/selscan -e ${DUMP}/logs/selscan \
    ${PIPE}/selscan.sh -w 100000 -g 200000 -s 20000 -x 1000000 \
      -h ${DUMP}/selection/ehh/${POP1}/${POP1}-${CHR}.hap \
      -r ${DUMP}/selection/ehh/${POP2}/${POP2}-${CHR}.hap \
      -m ${DUMP}/selection/ehh/CHR/${CHR}.map \
      -o ${DUMP}/selection/ehh/results/D-${POP1}-${POP2}-${CHR}
done

# XP-EHH post-processing
OUT=${DUMP}/selection/ehh/results
prefix="D-JFM-OTH-"
suffix=".xpehh.out"
NCBI=('NC_007070.3' 'NC_007071.3' 'NC_007072.3' \
'NC_007073.3' 'NC_007074.3' 'NC_007075.3' 'NC_007076.3' \
'NC_007077.3' 'NC_007078.3' 'NC_007079.3' 'NC_007080.3' \
'NC_007081.3' 'NC_007082.3' 'NC_007083.3' 'NC_007084.3' \
'NC_007085.3')
ENSEMBL=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16')

# Insert chromosome column at start of each file
x=0
for CHR in ${NCBI[@]}
do
  FILE=${OUT}/${prefix}${CHR}${suffix}
  cut -f2- ${FILE} > ${OUT}/test
#  sed -i '1s/^/chr\tpos\t/' ${OUT}/test
  sed "1 ! s/^/${ENSEMBL[${x}]}\t/" ${OUT}/test > ${FILE}
  x=$(($x+1))
done
# Create new field header
printf "chr\tstart\tgpos\tp1\tihh1\tp2\tihh2\txpehh\n" > ${OUT}/test
# Remove header from all out files
sed -i 1d ${OUT}/*.out
# Remove MT chromosomes
rm ${OUT}/*NC_001566*
# Merge all results together
cat ${OUT}/test ${OUT}/*.out > ${OUT}/${prefix}all${suffix}
# Clean up
rm test



# -----------------------------------------------------------------------------
# iHS
# -----------------------------------------------------------------------------

POP="JFM"
mkdir -p ${DUMP}/logs/selscan
for CHR in ${CHRS[@]}
do
  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs/selscan -e ${DUMP}/logs/selscan \
    ${PIPE}/ihs.sh \
      -h ${DUMP}/selection/ehh/${POP}/${POP}-${CHR}.hap \
      -m ${DUMP}/selection/ehh/CHR/${CHR}.map \
      -o ${DUMP}/selection/ehh/ihs/D-${POP}-${CHR}
done

# iHS post-processing
POP="OTH"
SUFFIX="out.100bins.norm"
ENSEMBL=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16')
x=0
for CHR in ${CHRS[@]}
do
  norm --files ${DUMP}/selection/ehh/ihs/D-${POP}-${CHR}.${SUFFIX}
  cp ${DUMP}/selection/ehh/ihs/D-${POP}-${CHR}.${SUFFIX} \
    ${DUMP}/selection/ehh/ihs/D-${POP}-${CHR}.${SUFFIX}.tmp
  sed "1 ! s/^/${ENSEMBL[${x}]}\t/" \
    ${DUMP}/selection/ehh/ihs/D-${POP}-${CHR}.${SUFFIX}.tmp \
    > ${DUMP}/selection/ehh/ihs/D-${POP}-${CHR}.${SUFFIX}
  rm ${DUMP}/selection/ehh/ihs/D-${POP}-${CHR}.${SUFFIX}.tmp
   x=$(($x+1))
done









#LG1    NC_007070.3
#LG2    NC_007071.3
#LG3    NC_007072.3
#LG4    NC_007073.3
#LG5    NC_007074.3
#LG6    NC_007075.3
#LG7    NC_007076.3
#LG8    NC_007077.3
#LG9    NC_007078.3
#LG10   NC_007079.3
#LG11   NC_007080.3
#LG12   NC_007081.3
#LG13   NC_007082.3
#LG14   NC_007083.3
#LG15   NC_007084.3
#LG16   NC_007085.3


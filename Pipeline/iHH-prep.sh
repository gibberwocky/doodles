#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# Load module
module load bioinfo/bcftools


# ==============================================================================
# iHH-prep.sh - a per-chromosome equivalent of - haps.sh
# ==============================================================================

# Clear variables 
OUT=
CHR=
FILE=
VCF=
POP=
START=
N=

while getopts ":p:c:o:v:i:x:n:" opt; do
  case $opt in
    p) OUT=${OPTARG};;
    c) CHR=${OPTARG};;
    o) FILE=${OPTARG};;
    v) VCF=${OPTARG};;
    i) POP=${OPTARG};;
    x) START=${OPTARG};;
    n) N=${OPTARG};;
  esac
done

if [[ -z ${OUT} ]] | [[ -z ${CHR} ]] | [[ -z ${FILE} ]] | [[ -z ${VCF} ]] | [[ -z ${POP} ]] | [[ -z ${START} ]] | [[ -z ${N} ]]
then
  usage
  exit 1
fi


# Extract chromosomes, excluding sites with uncalled genotypes, drop header
bcftools view -H -g ^het -g ^miss -O v --regions ${CHR} -o ${OUT}/${FILE}-${CHR} ${VCF}

# Extract header and save sample IDs
bcftools view -h  -O v -o ${OUT}/${CHR}-header.vcf ${VCF}
tail --l 1 ${OUT}/${CHR}-header.vcf  > ${OUT}/${CHR}-IDs
rm ${OUT}/${CHR}-header.vcf

# Prep directory and placeholder file
mkdir -p ${OUT}/${POP}
HAPFILE=${OUT}/${POP}/${POP}-${CHR}.hap
rm -f ${HAPFILE}
touch ${HAPFILE}

for((I=$((${START}));I<=$((${START}+${N}-1));I++));
do
  # Cut n paste
  cut -f${I} ${OUT}/${FILE}-${CHR} > ${OUT}/${CHR}-test.${I}
  cut -c 1 ${OUT}/${CHR}-test.${I} > ${OUT}/${CHR}-test.${I}a
  readarray tmp < ${OUT}/${CHR}-test.${I}a
  # Append data
  echo ${tmp[@]} >> ${HAPFILE}
  # Clean up
  rm ${OUT}/${CHR}-test.${I}
  rm ${OUT}/${CHR}-test.${I}a
done


# Extract mapfile details
mkdir -p ${OUT}/CHR
cut -f 1-3 ${OUT}/${FILE}-${CHR} | \
awk '{ temp = $1":"$2"[Amel4-5]"; $3 = temp; print $1"\t"$3"\t"$S2/28700"\t"$2 }' - > ${OUT}/CHR/${CHR}.map

# CHR cleanup
rm ${OUT}/${FILE}-${CHR}
rm ${OUT}/${CHR}-IDs


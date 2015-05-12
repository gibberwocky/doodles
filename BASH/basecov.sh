#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

REF=/save/dwragg/Apis/Apis_mellifera.fa
CHR=
OUT=
BAM=
Q=20

while getopts ":c:o:b:q:r:" opt; do
  case $opt in
    c) CHR=${OPTARG};;
    o) OUT=${OPTARG};;
    b) BAM=${OPTARG};;
    q) Q=${OPTARG};;
    r) REF=${OPTARG};;
  esac
done

if [[ -z ${CHR} ]] | [[ -z ${OUT} ]] | [[ -z ${BAM} ]]
then
  exit 1
fi

samtools mpileup -f ${REF} -q ${Q} -b ${BAM} -r ${CHR} -t DP -v -o ${OUT}


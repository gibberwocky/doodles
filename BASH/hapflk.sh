#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

IN=
OUT=
K=
CHR=

while getopts ":i:o:k:c::" opt; do
  case $opt in
    i) IN=${OPTARG};;
    o) OUT=${OPTARG};;
    k) K=${OPTARG};;
    c) CHR=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${OUT} ]] | [[ -z ${K} ]]         
then
  exit 1
fi


mkdir -p ${OUT}

# Run HAPFLK (on cluster!)
echo "Running HapFLK"
hapflk --bfile ${IN} \
  --chr ${CHR} \
  -K ${K} \
  --prefix ${OUT} \
  --phased





#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

PIPE=/save/seqapipop/Scripts
VEP=
VCF=
GTF=/save/dwragg/Apis/Apis_mellifera.GCA_000002195.1.25.gtf.clean
OUT=

while getopts ":e:v:g:o:" opt; do
  case $opt in
    e) VEP=${OPTARG};;
    v) VCF=${OPTARG};;
    g) GTF=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${VEP} ]] | [[ -z ${VCF} ]] | [[ -z ${OUT} ]]
then
  exit 1
fi

${PIPE}/cds.py -e ${VEP} -v ${VCF} -g ${GTF} -o ${OUT}


#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a


FAS=
OUT=

while getopts ":i:o:" opt; do
  case $opt in
    i) FAS=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${FAS} ]] | [[ -z ${OUT} ]]
then
  exit 1
fi

clustalw ${FAS} -outfile=${OUT} -output=nexus -quicktree
jmodeltest -d ${OUT} -w -g 11 -i -f -AIC -BIC -a -o ${OUT}.jmt


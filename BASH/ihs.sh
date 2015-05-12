#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

HAP=
MAP=
OUT=

while getopts ":h:m:o:" opt; do
  case $opt in
    h) HAP=${OPTARG};;
    m) MAP=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${HAP} ]] | [[ -z ${MAP} ]] | [[ -z ${OUT} ]]
then
  exit 1
fi

sed 's/2/1/g' ${HAP} > ${HAP}.01

selscan --ihs --trunc-ok --hap ${HAP}.01 --map ${MAP} --out ${OUT} 



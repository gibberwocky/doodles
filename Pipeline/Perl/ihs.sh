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

#selscan --ihs --skip-low-freq --trunc-ok --max-extend 100000 --hap ${HAP} --map ${MAP} --out ${OUT}

selscan --ihs --trunc-ok --ehh-win 20000 --hap ${HAP} --map ${MAP} --out ${OUT}



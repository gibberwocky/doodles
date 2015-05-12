#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

HAP=
REF=
MAP=
OUT=
MGAP=200000
SGAP=20000
WIN=100000
MAX=1000000

while getopts ":h:r:m:o:g:s:w:x:" opt; do
  case $opt in
    h) HAP=${OPTARG};;
    r) REF=${OPTARG};;
    m) MAP=${OPTARG};;
    o) OUT=${OPTARG};;
    g) MGAP=${OPTARG};;
    s) SGAP=${OPTARG};;
    w) WIN=${OPTARG};;
    x) MAX=${OPTARG};;
  esac
done

if [[ -z ${HAP} ]] | [[ -z ${REF} ]] | [[ -z ${MAP} ]] | [[ -z ${OUT} ]]
then
  exit 1
fi

selscan --xpehh --trunc-ok --ehh-win ${WIN} --max-gap ${MGAP} --max-extend ${MAX} --gap-scale ${SGAP} --hap ${HAP} --ref ${REF} --map ${MAP} --out ${OUT}

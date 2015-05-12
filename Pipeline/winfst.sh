#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# selection.sh
# ==============================================================================


# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
#RSCRIPTS=/home/dwragg/save/Scripts/R
RSCRIPTS=/home/dwragg/work/Pipeline
CHRINFO=/save/dwragg/Apis/chrInfo
IN=
PED=
MAP=
PHENO=
CHR=
WIN=
OLAP=1000

while getopts ":i:p:m:h:c:w:f:o:x:" opt; do
  case $opt in
    i) IN=${OPTARG};;
    p) PED=${OPTARG};;
    m) MAP=${OPTARG};;
    h) PHENO=${OPTARG};;
    c) CHR=${OPTARG};;
    w) WIN=${OPTARG};;
    f) CHRINFO=${OPTARG};;
    o) OUT=${OPTARG};;
    x) OLAP=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${PED} ]] | [[ -z ${MAP} ]] | [[ -z ${PHENO} ]] | [[ -z ${CHR} ]] | [[ -z ${WIN} ]]
then
  usage
  exit 1
fi

${RSCRIPTS}/winfst.R ${IN} ${PED} ${MAP} ${PHENO} ${CHR} ${WIN} ${CHRINFO} ${OUT} ${OLAP}

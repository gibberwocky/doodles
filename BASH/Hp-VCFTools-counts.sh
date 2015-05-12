#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# Hp-VCFtools-counts.sh
# ==============================================================================


# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
RSCRIPTS=/home/dwragg/save/Scripts/R
CHRINFO=/save/dwragg/Apis/chrInfo
INDIR=
IN=
WIN=5000
OVERLAP=1000

while getopts ":d:i:w:x:" opt; do
  case $opt in
    d) INDIR=${OPTARG};;
    i) IN=${OPTARG};;
    w) WIN=${OPTARG};;
    x) OVERLAP=${OPTARG};;
  esac
done

if [[ -z ${INDIR} ]] | [[ -z ${IN} ]] 
then
  usage
  exit 1
fi

${RSCRIPTS}/Hp-VCFTools-counts.R ${INDIR} ${IN} ${WIN} ${OVERLAP} ${CHRINFO}

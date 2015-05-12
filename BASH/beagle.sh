#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a


# ==============================================================================
# BEAGLE
# ==============================================================================

IN=
OUT=
BEAGLE="/usr/local/bioinfo/src/Beagle/beagle4.jar"

while getopts ":i:o:" opt; do
  case $opt in
    i) IN=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${OUT} ]]
then
  exit 1
fi


# Phase and impute missing data
java -d64 -jar ${BEAGLE} gt=${IN} out=${OUT} usephase=true


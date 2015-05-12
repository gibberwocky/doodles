#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

OUT=
IDFILE=
PHASE=

while getopts ":o:i:p:" opt; do
  case $opt in
    o) OUT=${OPTARG};;
    i) IDFILE=${OPTARG};;
    p) PHASE=${OPTARG};;
  esac
done

if [[ -z ${OUT} ]] | [[ -z ${IDFILE} ]] | [[ -z ${PHASE} ]]
then
  exit 1
fi

fs ${OUT}.cp -idfile ${IDFILE} -phasefiles ${PHASE} -go

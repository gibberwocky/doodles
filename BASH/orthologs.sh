#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a


# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

IN=
DB=
FORMAT=8
TOOL="blastx"
EVAL=0.5
MEGA="T"
OUT=

while getopts ":i:d:f:t:e:m:o:" opt; do
  case $opt in
    i) IN=${OPTARG};;
    d) DB=${OPTARG};;
    f) FORMAT=${OPTARG};;
    t) TOOL=${OPTARG};;
    e) EVAL=${OPTARG};;
    m) MEGA=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${DB} ]] | [[ -z ${OUT} ]]
then
  usage
  exit 1
fi

# BLAST!
blastall -i ${IN} -d ${DB} -m ${FORMAT} -p ${TOOL} -e ${EVAL} -n ${MEGA} > ${OUT}


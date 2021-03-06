#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# vep.sh
# ==============================================================================


# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
IN=
OUT=

while getopts ":i:o:" opt; do
  case $opt in
    i) IN=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${OUT} ]]  
then
  usage
  exit 1
fi

# 76 registry
# /usr/local/bioinfo/src/ensembl-api/branch-76/variant_effect_predictor.registry

# 77 registry
# /usr/local/bioinfo/src/ensembl-api/current/variant_effect_predictor.registry

variant_effect_predictor.pl \
  --database \
  --species bee \
  --registry /usr/local/bioinfo/src/ensembl-api/branch-78/variant_effect_predictor.registry \
  --input_file ${IN} \
  --output_file ${OUT}_vep_output.txt \
  --symbol \
  --stats_text \
  --force_overwrite



#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

BEAGLE=/usr/local/bioinfo/src/Beagle/beagle4.jar
VCF=
OUT=

while getopts ":v:o:" opt; do
  case $opt in
    v) VCF=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${VCF} ]] | [[ -z ${OUT} ]] 
then
  exit 1
fi

java -d64 -jar ${BEAGLE} gt=${VCF} out=${OUT}
tabix -fp vcf ${OUT}.vcf.gz


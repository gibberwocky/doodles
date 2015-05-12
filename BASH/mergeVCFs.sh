#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

PIPE=/save/seqapipop/Scripts
POP=
IN=
OUT=

while getopts ":p:i:o:" opt; do
  case $opt in
    p) POP=${OPTARG};;
    i) IN=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${POP} ]] | [[ -z ${OUT} ]] | [[ -z ${IN} ]]
then
  exit 1
fi




# Merge diploids into single VCF
bcftools merge -m snps -O z \
  -l ${IN}/${POP}-diploid_vcfs.list \
  -o ${OUT}/${POP}_diploid.vcf.gz
tabix -fp vcf ${OUT}/${POP}_diploid.vcf.gz

# Remove unnecessary files
readarray FILES < ${IN}/${POP}-diploid_vcfs.list
rm ${FILES[@]}
tmp=$(printf "%s.tbi " ${FILES[@]})
rm ${tmp}



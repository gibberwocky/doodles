#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

PIPE=/save/seqapipop/Scripts
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

cd ${OUT}
bcftools view --type snps --exclude 'AN!=2' -O z -o ${VCF}-clean.vcf.gz ${VCF}.vcf.gz
bcftools annotate -x FORMAT -O z -o ${VCF}.vcf.gz ${VCF}-clean.vcf.gz
tabix -fp vcf ${VCF}.vcf.gz
rm ${VCF}-clean.vcf.gz



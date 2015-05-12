#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

BAYESCAN=" /usr/local/bioinfo/src/BayeScan/current/binaries"
PGDSpider="/usr/local/bioinfo/src/PGDSspider/current"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
IN=
OUT=
CHR=1
SPID=
FAM=

while getopts ":i:o:c:s:f:" opt; do
  case $opt in
    i) IN=${OPTARG};;
    o) OUT=${OPTARG};;
    c) CHR=${OPTARG};;
    s) SPID=${OPTARG};;
    f) FAM=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${OUT} ]] | [[ -z ${CHR} ]]  | [[ -z ${SPID} ]] | [[ -z ${FAM} ]]        
then
  exit 1
fi


mkdir -p ${OUT}/in
mkdir -p ${OUT}/out
mkdir -p ${OUT}/logs

i=1
#for((i=1;i<=${CHR};i++))
#do

  # First output to PED format
#  echo "Writing PED file: Chr ${i}"
#  ${PLINK}/plink --bfile ${IN} \
#    --fam ${FAM} \
#    --allow-no-sex \
#    --allow-extra-chr \
#    --chr-set 17 \
#    --chr ${i} \
#    --set-missing-snp-ids @:#[Amel4-5] \
#    --geno 0.5 \
#    --recode 12 \
#    --bp-space 100 \
#    --out ${OUT}/in/${i}

  # Convert PED to BAYESCAN format
#  echo "Converting PED to BayeScan file: Chr ${i}"
#  java -jar ${PGDSpider}/PGDSpider2-cli.jar \
#    -inputfile ${OUT}/in/${i}.ped \
#    -inputformat PED \
#    -outputfile ${OUT}/in/${i}.bayes \
#    -outputformat GESTE_BAYE_SCAN \
#    -spid ${SPID}

  # Run BAYESCAN (on cluster!)
  echo "Running BayeScan: Chr ${i}"
  ${BAYESCAN}/BayeScan2.1_linux64bits \
    ${OUT}/in/${i}.bayes \
    -od ${OUT}/out \
    -pr_idds 1000 \
    -threads 16

#done





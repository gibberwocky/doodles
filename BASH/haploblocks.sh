#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# admixture.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Block analyses\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -b <plink bfile> -c <chromosome> -o <output path and name> [options]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script runs Plink's haploblock analysis.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-f <chr>\tPath to Plink FAM file\e[0m"
  echo -e "\e[92m\t-m <float>\tMinor allele frequency (MAF) threshold\e[0m"
  echo -e "\e[92m\t-k <int>\tMaximum Kb for block calculation\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mPlink\e[0m"
  echo -e "\n"

}

PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
BFILE=
CHR=
OUT=
FAM=
MAF=0.05
MAXKB=1000

while getopts ":b:c:o:f:m:k:" opt; do
  case $opt in
    b) BFILE=${OPTARG};;
    c) CHR=${OPTARG};;
    o) OUT=${OPTARG};;
    f) FAM=${OPTARG};;
    m) MAF=${OPTARG};;
    k) MAXKB=${OPTARG};;
  esac
done

if [[ -z ${BFILE} ]] | [[ -z ${CHR} ]] | [[ -z ${OUT} ]]
then
  usage
  exit 1
fi

${PLINK}/plink --bfile ${BFILE} \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr ${CHR} \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${OUT} \
  --blocks no-pheno-req \
  --blocks-min-maf ${MAF} \
  --blocks-max-kb ${MAXKB} \
  --keep ${FAM}






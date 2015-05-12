#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# Set-up
module load bioinfo/bcftools

# ==============================================================================
# haps2.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Haplotype alleles\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -v <VCF file> -o <output path> -p <population name> -s <sample names file> -c <chromosome>\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script creates hap and map files for a population within a VCF file.\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mBCFtools\e[0m"
  echo -e "\n"

}


VCFT="/usr/local/bioinfo/src/vcftools/vcftools_0.1.12a/bin"
VCF=
OUT=
POP=
SAMPLES=
CHR=

while getopts ":v:o:p:s:c:" opt; do
  case $opt in
    v) VCF=${OPTARG};;
    o) OUT=${OPTARG};;
    p) POP=${OPTARG};;
    s) SAMPLES=${OPTARG};;
    c) CHR=${OPTARG};;
  esac
done

if [[ -z ${VCF} ]] | [[ -z ${OUT} ]] | [[ -z ${POP} ]] | [[ -z ${SAMPLES} ]] | [[ -z ${CHR} ]]
then
  usage
  exit 1
fi

# Convert VCF to haplotype format
${VCFT}/vcftools --gzvcf ${VCF} --012 --chr ${CHR} --keep ${SAMPLES} --out ${OUT}/${POP}-${CHR}

# Substitute tabs for spaces
sed -i 's/\t/ /g' ${OUT}/${POP}-${CHR}.012

# Substitute 2's for 1's
#sed -i 's/2/1/g' ${OUT}/${POP}-${CHR}.012

# Remove first column (sample ID) and move to population folder
cut -d" " -f2- ${OUT}/${POP}-${CHR}.012 > ${OUT}/${POP}/${POP}-${CHR}.hap

# Move to own folder
mv ${OUT}/${POP}-${CHR}.012.indv ${OUT}/${POP}/${POP}-${CHR}.ind

# Create map file
if [ -e "${OUT}/CHR/${CHR}.map" ]
then
  echo "${OUT}/CHR/${CHR}.map exists, not re-writing"
else
  awk '{ print $1"\t"$1":"$2"\t"$2/27027"\t"$2 }' ${OUT}/${POP}-${CHR}.012.pos > ${OUT}/CHR/${CHR}.map
fi

# Clean up
rm ${OUT}/${POP}-${CHR}.012*



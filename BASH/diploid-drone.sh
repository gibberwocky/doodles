#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a
module load bioinfo/bcftools

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e


# ==============================================================================
# diploid-drone.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Create faux-diploid\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -v <VCF file containing haploids> -a <sample 1 ID> -b <sample 2 ID> -o <output path> -s <output name>\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script generates faux-diploids from haploid SNP data.\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mdiploid-drone-v2.pl\e[0m"
  echo -e "\e[92mbcftools\e[0m"
  echo -e "\n"

}



PERL="/save/seqapipop/Scripts/Perl"
V=
A=
B=
OUT=
ID=

while getopts ":v:a:b:s:o:" opt; do
  case $opt in
    v) V=${OPTARG};;
    a) A=${OPTARG};;
    b) B=${OPTARG};;
    o) OUT=${OPTARG};;
    s) ID=${OPTARG};;
  esac
done

if [[ -z ${V} ]] | [[ -z ${A} ]] | [[ -z ${B} ]] | [[ -z ${OUT} ]] | [[ -z ${ID} ]]
then
  usage
  exit 1
fi



# Extract a pair of drones to vcf
bcftools view --types snps -M 2 -O v \
  -s ${A},${B} \
  -o ${OUT}/${ID} \
  ${V}

# Create diploid drone
perl ${PERL}/diploid-drone-v3.pl ${OUT}/${ID} ${ID} > ${OUT}/${ID}-fix.vcf
rm ${OUT}/${ID}
sed -i '/^$/d' ${OUT}/${ID}-fix.vcf

# Extract drone from VCF
bcftools view --types snps -M 2 -O v \
  -s ${ID} \
  -o ${OUT}/${ID}.vcf \
  ${OUT}/${ID}-fix.vcf
rm ${OUT}/${ID}-fix.vcf

# Compress and index
bgzip -f ${OUT}/${ID}.vcf
tabix -p vcf ${OUT}/${ID}.vcf.gz


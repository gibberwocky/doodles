#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# Modify environment for Java7
module load bioinfo/Java7 

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# mergeBAMs-v2.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Merge multiple BAMs\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -o <output path> -n <output name> -b <file containing list of BAMs to merge>\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script merges multiple BAM files.\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mSamtools\e[0m"
  echo -e "\n"

}

# Variables
NAME=
OUT=
BAMS=

while getopts ":n:b:o:" opt; do
  case $opt in
    n) NAME=${OPTARG};;
    b) BAMS=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${NAME} ]] | [[ -z ${BAMS} ]] | [[ -z ${OUT} ]]
then
  usage
  exit 1
fi


# ==============================================================================
# Merge BAM files
# ==============================================================================
ID=($(cut -f1 ${BAMS}))
samtools merge ${OUT}/${NAME}.tmp.bam ${ID[@]}
samtools sort ${OUT}/${NAME}.tmp.bam ${OUT}/${NAME}
rm ${OUT}/${NAME}.tmp.bam
samtools index ${OUT}/${NAME}.bam



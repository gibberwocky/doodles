#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# downsample.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Downsample\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -i <input BAM file> -o <output BAM file> -p <probability of read being retained>\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script downsamples a BAM file by retaining reads at a given probability.\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mPicard\e[0m"
  echo -e "\n"

}

# Clear variables 
INFILE=
PROB=
OUT=
PICARD="/usr/local/bioinfo/src/picard-tools/current"

while getopts ":i:o:p:" opt; do
  case $opt in
    i) INFILE=${OPTARG};;
    o) OUT=${OPTARG};;
    p) PROB=${OPTARG};;
  esac
done

if [[ -z ${INFILE} ]] | [[ -z ${OUT} ]] | [[ -z ${PROB} ]]
then
  usage
  exit 1
fi


# Downsample BAM/SAM
java -d64 -jar ${PICARD}/DownsampleSam.jar \
  INPUT=${INFILE} \
  OUTPUT=${OUT} \
  PROBABILITY=${PROB} \
  QUIET=T \
  VERBOSITY=ERROR \
  VALIDATION_STRINGENCY=LENIENT


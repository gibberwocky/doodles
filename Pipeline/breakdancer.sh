#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# Load Java (just in case)
module load bioinfo/Java7


# ==============================================================================
# breakdancer.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: BreakDancer analyses\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -a <path to BAM file 1> -b <path to BAM file 2> -o <output path> -s <output name> [options]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script runs BreakDancer.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-q <int>\tQuality threshold [default = 20]\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mBreakDancer\e[0m"
  echo -e "\n"

}

BDANCE="/usr/local/bioinfo/src/BreakDancer/breakdancer-1.3.6/bin"
BAM2CFG="/usr/local/bioinfo/src/BreakDancer/breakdancer-1.3.6/perl"
BAM1=
BAM2=
OUT=
NAME=
Q=20

while getopts ":a:b:o:s:q:" opt; do
  case $opt in
    a) BAM1=${OPTARG};;
    b) BAM2=${OPTARG};;
    o) OUT=${OPTARG};;
    s) NAME=${OPTARG};;
    q) Q=${OPTARG};;
  esac
done

if [[ -z ${BAM1} ]] | [[ -z ${BAM2} ]] | [[ -z ${OUT} ]] | [[ -z ${NAME} ]]
then
  usage
  exit 1
fi


cd ${OUT}

echo -e "\e[91mCreating config file from the BAM files\e[94m"
#perl ${BAM2CFG}/bam2cfg.pl -q ${Q} -g -h ${BAM1} ${BAM2} > ${OUT}/${NAME}.cfg

echo -e "\e[91mCalculationg structural variants\e[94m"
mkdir -p ${OUT}/fastq
${BDANCE}/breakdancer-max -fh -q ${Q} -d ${OUT}/fastq/${NAME} -g ${OUT}/${NAME}.bed ${OUT}/${NAME}.cfg > ${OUT}/${NAME}.ctx




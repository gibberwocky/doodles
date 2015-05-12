#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# Modify environment for Java7
module load bioinfo/Java7 

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# stat.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: GC content by interval\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
}

# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
REF=/home/dwragg/save/Apis/Apis_mellifera.fa
IN=
OUT=
INTERVALS=

while getopts ":i:l:o:" opt; do
  case $opt in
    i) IN=${OPTARG};;
    l) INTERVALS=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${INTERVALS} ]] | [[ -z ${OUT} ]]
then
  usage
  exit 1
fi





# Depth of coverage via GATK
echo -e "\e[91m[GATK:GCContentByInterval \n\t\e[94m\e[1m<- \e[0m\e[92m/${IN} \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${IN}_GC\e[0m"
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T GCContentByInterval \
  -R ${REF} \
  -I ${IN} \
  -o ${OUT}/${IN}_GC \
  -L ${INTERVALS} \
  -l FATAL 














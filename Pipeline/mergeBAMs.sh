#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# Modify environment for Java7
module load bioinfo/Java7 

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# mergeBAMs.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m==============================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Merge and pre-process two BAM files\e[0m"
echo -e "\e[1m\e[34m==============================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -s <sample> -a <path and name of BAM file 1> -b <path and name of BAM file 2> -o <output path> [options]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script merges two BAM files, replaces RG tags, sorts, indexes, and marks duplicates.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-i <str>\tPath to input file containing pipeline parameters [defafult = /save/seqapipop/Scripts/params]\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mPicard\e[0m"
  echo -e "\e[92mSamtools\e[0m"
  echo -e "\n"

}

# Variables
INFILE=/save/seqapipop/Scripts/params
IN=
OUT=
RUN1=
RUN2=

while getopts ":i:s:a:b:o:" opt; do
  case $opt in
    i) INFILE=${OPTARG};;
    s) IN=${OPTARG};;
    a) RUN1=${OPTARG};;
    b) RUN2=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${INFILE} ]] | [[ -z ${IN} ]] | [[ -z ${RUN1} ]] | [[ -z ${RUN2} ]] | [[ -z ${OUT} ]]
then
  usage
  exit 1
fi

while IFS="=" read name value
do

  if [[ "$name" == "REF" ]]; then REF=${value//\"/}; fi
  if [[ "$name" == "SAMPLE" ]]; then SAMPLE=${value//\"/}; fi
  if [[ "$name" == "BWA" ]]; then BWA=${value//\"/}; fi
  if [[ "$name" == "PICARD" ]]; then PICARD=${value//\"/}; fi
  if [[ "$name" == "GATK" ]]; then GATK=${value//\"/}; fi
  if [[ "$name" == "PLATYPUS" ]]; then PLATYPUS=${value//\"/}; fi
  if [[ "$name" == "VCFT" ]]; then VCFT=${value//\"/}; fi
  if [[ "$name" == "BAYSIC" ]]; then BAYSIC=${value//\"/}; fi
  if [[ "$name" == "PERLSCRIPTS" ]]; then PERLSCRIPTS=${value//\"/}; fi
  if [[ "$name" == "RSCRIPTS" ]]; then RSCRIPTS=${value//\"/}; fi
  if [[ "$name" == "BEAGLE" ]]; then BEAGLE=${value//\"/}; fi

done < ${INFILE}



# ==============================================================================
# Merge BAM files
# ==============================================================================
java -d64 -jar $PICARD/MergeSamFiles.jar \
  INPUT=${RUN1} \
  INPUT=${RUN2} \
  OUTPUT=${OUT}/${IN}/${IN}_tmp.bam \
  SORT_ORDER=coordinate \
  USE_THREADING=T \
  QUIET=T \
  VERBOSITY=ERROR \
  VALIDATION_STRINGENCY=LENIENT

# Index BAM file
samtools index ${OUT}/${IN}/${IN}_tmp.bam

# Replace RG information
java -d64 -jar $PICARD/AddOrReplaceReadGroups.jar \
  INPUT=${OUT}/${IN}/${IN}_tmp.bam \
  OUTPUT=${OUT}/${IN}/${IN}_bootstrap.bam \
  SORT_ORDER=coordinate \
  RGID="SeqApiPop" \
  RGLB=${IN} \
  RGPL="illumina" \
  RGPU="INRA" \
  RGSM=${IN} \
  QUIET=T \
  VERBOSITY=ERROR \
  VALIDATION_STRINGENCY=LENIENT

# Remove tmp files
rm ${OUT}/${IN}/${IN}_tmp.*

# Index BAM file
samtools index ${OUT}/${IN}/${IN}_bootstrap.bam

# Mark Duplicates
java -d64 -jar ${PICARD}/MarkDuplicates.jar \
  INPUT=${OUT}/${IN}/${IN}_bootstrap.bam \
  OUTPUT=${OUT}/${IN}/${IN}_bootdup.bam \
  METRICS_FILE=${OUT}/${IN}/metrics/${IN}_dup.metrics \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
  QUIET=T \
  VERBOSITY=ERROR \
  VALIDATION_STRINGENCY=LENIENT

# Remove dup file
rm ${OUT}/${IN}/${IN}_bootdup.bam


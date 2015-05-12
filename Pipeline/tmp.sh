#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# Load environment
module load bioinfo/Java7

REF=/home/dwragg/save/Apis/Apis_mellifera.fa
GATK=/usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.1-1
DUMP="/home/dwragg/work/Analysis"

java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${REF} \
  -I ${DUMP}/vcfs/lists/Harpur-SeqApiPop_bams.list \
  -o ${DUMP}/vcfs/GATKcov \
  -ct 2 -ct 5 -ct 8 \
  --omitDepthOutputAtEachBase \
  --omitIntervalStatistics \
  --omitLocusTable \
  -l FATAL 


#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a


# ==============================================================================
# ANGSD notes
# ==============================================================================

# Working directory
DUMP="/home/dwragg/work/tmp"

# Output list of JFM population BAM files to file (JFM.bamlist)
JFM=(${DUMP}/JFM*/*.bam)
ls ${JFM[@]} > ${DUMP}/angsd/JFM.bamlist

# Phase and impute missing data
BEAGLE="/usr/local/bioinfo/src/Beagle"
java -d64 -jar ${BEAGLE}/beagle3.jar \
  like=${DUMP}/angsd/JFM.beagle.gz \
  out=${DUMP}/beagle/JFM


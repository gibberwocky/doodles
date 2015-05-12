#!/bin/bash
# IEKNGHER1ujm

ln -fs /work/vignal/Analysis/Project_ROYALBEE.301/Run_Lot2_Pool1_align.6437/RawData ${DUMP}



# Set root paths
PIPE="/save/seqapipop/Scripts"

DUMP="/home/dwragg/work"
IN=/work/dwragg/Project_ROYALBEE.301/Run_Pool7.6708/RawData
ID="AOC32_GATCAG_L003"
ID="AOC33_TAGCTT_L003"

# file name format of fastq files: $ID_R1.fastq.gz $ID_R2.fastq.gz

# map.sh	stat.sh		snps.sh
  mkdir -p ${DUMP}/logs/${ID}
  cd ${DUMP}
  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs/${ID} \
    -e ${DUMP}/logs/${ID} \
    ${PIPE}/snps.sh -s ${ID} -f ${IN} -o ${DUMP}/SRA






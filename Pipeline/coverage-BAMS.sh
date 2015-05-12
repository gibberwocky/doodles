# Modify environment for Java7
module load bioinfo/Java7 

REF=/home/dwragg/save/Apis/Apis_mellifera.fa
GATK=/usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.1-1
DUMP="/home/dwragg/work/Analysis"

# BAM arrays
readarray BAMS_A < ${DUMP}/vcfs/lists/A_bams.list
readarray BAMS_M < ${DUMP}/vcfs/lists/M_bams.list
readarray BAMS_C < ${DUMP}/vcfs/lists/C_bams.list
readarray BAMS_O < ${DUMP}/vcfs/lists/O_bams.list
readarray BAMS_JFM < ${DUMP}/vcfs/lists/JFM_bams.list
readarray BAMS_OTH < ${DUMP}/vcfs/lists/OTH_bams.list
tmp=(${BAMS_A[@]} ${BAMS_M[@]} ${BAMS_C[@]} ${BAMS_O[@]} ${BAMS_JFM[@]} ${BAMS_OTH[@]})

printf "%s\n" "${tmp[@]}" > ${DUMP}/vcfs/lists/Harpur-SeqApiPop_bams.list

# GATK calculate depth of coverage
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



# Above script hashed into tmp.sh to run on cluster as is time-consuming
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/tmp.sh

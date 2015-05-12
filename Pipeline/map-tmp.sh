#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# Modify environment for Java7
module load bioinfo/Java7 

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# MAP.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m======================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop \e[0m\e[93mNGS pipeline\e[0m"
echo -e "\e[1m\e[34m======================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -i <parameters file> -s <sample> [options]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis pipeline has been developed to map Illumina HiSeq reads generated for the SeqApiPop project. All of the parameters in the parameter file require setting. This module requires the installation of BWA, GATK, and Picard.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-i <str>\tPath to input file containing pipeline parameters\e[0m"
  echo -e "\e[92m\t-f <str>\tPath to folder containing fastq.gz files to be mapped\e[0m"
  echo -e "\e[92m\t-s <str>\tSample read file name prefix\e[0m"
  echo -e "\e[92m\t-b <int>\tNumber of bootstrap iterations during BQSR [3]\e[0m"
  echo -e "\n"

}

# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
INFILE=
PAIRED="T"
QFIX="F"
BOOTS=2
REF=
IN=
OUT=
SAMPLE=
BWA=
PICARD=
GATK=
PLATYPUS=
VCFT=
PERLSCRIPTS=
RSCRIPTS=
BAYSIC=

while getopts ":i:p:f:s:b:o:q:" opt; do
  case $opt in
    i) INFILE=${OPTARG};;
    f) IN=${OPTARG};;
    s) SAMPLE=${OPTARG};;
    b) BOOTS=${OPTARG};;
    o) OUT=${OPTARG};;
    p) PAIRED=${OPTARG};;
    q) QFIX=${OPTARG};;
  esac
done

if [[ -z ${INFILE} ]] | [[ -z ${SAMPLE} ]] | [[ -z ${IN} ]] | [[ -z ${OUT} ]]
then
  usage
  exit 1
fi

echo -e "\n\e[91mScript variables:\e[0m"
echo -e "\e[96mParameters file: \e[92m${INFILE}\e[0m"
echo -e "\e[96mBQSR iterations: \e[92m${BOOTS}\e[0m"


# ==============================================================================
# Read in variables
# ==============================================================================

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

done < ${INFILE}

if [[ "$REF" == "" ]]; then echo "REF missing"; exit 0; fi
if [[ "$OUT" == "" ]]; then echo "OUT path missing"; exit 0; fi
if [[ "$SAMPLE" == "" ]]; then echo "SAMPLE missing"; exit 0; fi
if [[ "$BWA" == "" ]]; then echo "BWA missing"; exit 0; fi
if [[ "$PICARD" == "" ]]; then echo "PICARD path missing"; exit 0; fi
if [[ "$GATK" == "" ]]; then echo "GATK path missing"; exit 0; fi
if [[ "$PLATYPUS" == "" ]]; then echo "PLATYPUS path missing"; exit 0; fi
if [[ "$VCFT" == "" ]]; then echo "VCFT path missing"; exit 0; fi
if [[ "$BAYSIC" == "" ]]; then echo "BAYSIC path missing"; exit 0; fi
if [[ "$PERLSCRIPTS" == "" ]]; then echo "PERLSCRIPTS path missing"; exit 0; fi
if [[ "$RSCRIPTS" == "" ]]; then echo "RSCRIPTS path missing"; exit 0; fi

echo -e "\n\e[91mVariables imported succesfully:\e[0m"

echo -e "\n\e[96mSAMPLE: \e[92m${SAMPLE}\e[0m"
echo -e "\e[96mPAIRED: \e[92m${PAIRED}\e[0m"
echo -e "\e[96mREF: \e[92m${REF}\e[0m"
echo -e "\e[96mIN: \e[92m${IN}\e[0m"
echo -e "\e[96mOUT: \e[92m${OUT}\e[0m"
echo -e "\e[96mBWA: \e[92m${BWA}\e[0m"
echo -e "\e[96mPICARD: \e[92m${PICARD}\e[0m"
echo -e "\e[96mGATK: \e[92m${GATK}\e[0m"
echo -e "\e[96mPLATYPUS: \e[92m${PLATYPUS}\e[0m"
echo -e "\e[96mVCFT: \e[92m${VCFT}\e[0m"
echo -e "\e[96mBAYSIC: \e[92m${BAYSIC}\e[0m"
echo -e "\e[96mPERLSCRIPTS: \e[92m${PERLSCRIPTS}\e[0m"
echo -e "\e[96mRSCRIPTS: \e[92m${RSCRIPTS}\e[0m"


# Create folders for storing files
echo -e "\e[91m\nCreating folders:\e[0m"
echo -e "\e[92m	${OUT}/${SAMPLE}/logs\e[0m"
mkdir -p ${OUT}/${SAMPLE}/logs
echo -e "\e[92m	${OUT}/${SAMPLE}/metrics\e[0m"
mkdir -p ${OUT}/${SAMPLE}/metrics
echo -e "\e[92m	${OUT}/${SAMPLE}/vcfs\e[0m"
mkdir -p ${OUT}/${SAMPLE}/vcfs


# Number of steps in pipeline, this will need updating if additional steps are added, to provide an idea of how long is remaining whilst running the pipe
STEPS=7



# ==============================================================================
# Mapping
# ==============================================================================

# Map Reads
#if [ "${PAIRED}" = "T" ]; then
#  echo -e "\e[91m[1/${STEPS}] BWA:Map reads [BWA MEM] \n\t\e[94m\e[1m<- \e[0m\e[92m${REF} \n\t\e[94m\e[1m<- \e[0m\e[92m${IN}/${READS_1} \n\t\e[94m\e[1m<- \e[0m\e[92m${IN}/${READS_2} \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}.sam \e[0m"
#  READS_1=${SAMPLE}_R1.fastq.gz
#  READS_2=${SAMPLE}_R2.fastq.gz
##  ${BWA}/bwa mem -M -R @RG"\t"ID:SeqApiPop"\t"SM:${SAMPLE}"\t"PL:illumina"\t"LB:${SAMPLE}"\t"PU:INRA ${REF} ${IN}/${READS_1} ${IN}/${READS_2} > ${OUT}/${SAMPLE}/${SAMPLE}.sam 
#fi
#if [ "${PAIRED}" = "F" ]; then
#  echo -e "\e[91m[1/${STEPS}] BWA:Map reads [BWA MEM] \n\t\e[94m\e[1m<- \e[0m\e[92m${REF} \n\t\e[94m\e[1m<- \e[0m\e[92m${IN}/${READS_1} \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}.sam \e[0m"
#  READS_1=${SAMPLE}.fastq.gz
##  ${BWA}/bwa mem -M -R @RG"\t"ID:SeqApiPop"\t"SM:${SAMPLE}"\t"PL:illumina"\t"LB:${SAMPLE}"\t"PU:INRA ${REF} ${IN}/${READS_1} > ${OUT}/${SAMPLE}/${SAMPLE}.sam 
#fi

# Sort SAM into coordinate order and save as BAM
#echo -e "\e[91m[2/${STEPS}] Picard:SortSam [coordinate] \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}.sam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}.bam\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardSortSam.log
#java -d64 -jar $PICARD/SortSam.jar \
#  INPUT=${OUT}/${SAMPLE}/${SAMPLE}.sam \
#  OUTPUT=${OUT}/${SAMPLE}/${SAMPLE}.bam \
#  SORT_ORDER=coordinate \
#  QUIET=T \
#  VERBOSITY=ERROR \
#  VALIDATION_STRINGENCY=LENIENT \
#  2> >(tee "$logfile")
#rm ${OUT}/${SAMPLE}/${SAMPLE}.sam

# Mark duplicates 
#echo -e "\e[91m[3/${STEPS}] Picard:MarkDuplicates \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}_dup.bam\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardMarkDuplicates.log
#java -d64 -jar ${PICARD}/MarkDuplicates.jar \
#  INPUT=${OUT}/${SAMPLE}/${SAMPLE}.bam \
#  OUTPUT=${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
#  METRICS_FILE=${OUT}/${SAMPLE}/metrics/${SAMPLE}_dup.metrics \
#  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
#  QUIET=T \
#  VERBOSITY=ERROR \
#  VALIDATION_STRINGENCY=LENIENT \
#  2> >(tee "$logfile")
#rm ${OUT}/${SAMPLE}/${SAMPLE}.bam

# Index BAM file 
#echo -e "\e[91m[4/${STEPS}]:BuildBamIndex\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardBuildBamIndex.log
#java -d64 -jar ${PICARD}/BuildBamIndex.jar \
#  INPUT=${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
#  QUIET=T \
#  VERBOSITY=ERROR \
#  2> >(tee "$logfile")



# ==============================================================================
# Local realignment around INDELs
# ==============================================================================

# Create target interval list
#echo -e "\e[91m[5/${STEPS}] GATK:RealignerTargetCreator \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/logs/${SAMPLE}_dup_intervals.list\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_RealignerTargetCreator.err

#if [ "${QFIX}" = "T" ]; then
#  java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#    -T RealignerTargetCreator \
#    -R ${REF} \
#    -I ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
#    -o ${OUT}/${SAMPLE}/logs/${SAMPLE}_dup_intervals.list \
#    -l FATAL \
#    --allow_potentially_misencoded_quality_scores \
#    2> >(tee "$logfile")
#fi

#if [ "${QFIX}" = "F" ]; then
#  java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#    -T RealignerTargetCreator \
#    -R ${REF} \
#    -I ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
#    -o ${OUT}/${SAMPLE}/logs/${SAMPLE}_dup_intervals.list \
#    -l FATAL \
#    2> >(tee "$logfile")
#fi

# Perform realignment around INDELs 
#echo -e "\e[91m[6/${STEPS}] GATK:IndelRealigner \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_IndelRealigner.err

#if [ "${QFIX}" = "T" ]; then
#java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#  -T IndelRealigner \
#  -R ${REF} \
#  -I ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
#  -targetIntervals ${OUT}/${SAMPLE}/logs/${SAMPLE}_dup_intervals.list \
#  -o ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
#  -l FATAL \
#  --allow_potentially_misencoded_quality_scores \
#  2> >(tee "$logfile")
#fi

#if [ "${QFIX}" = "F" ]; then
#java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#  -T IndelRealigner \
#  -R ${REF} \
#  -I ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
#  -targetIntervals ${OUT}/${SAMPLE}/logs/${SAMPLE}_dup_intervals.list \
#  -o ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
#  -l FATAL \
#  2> >(tee "$logfile")
#fi


# Delete surplus files
#rm ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam
#rm ${OUT}/${SAMPLE}/${SAMPLE}_dup.bai



# ==============================================================================
# BQSR based on boostrapping UnifiedGenotyper and BaseRecalibrator
# ==============================================================================

#echo -e "\e[91m[7/${STEPS}] GATK:Base Quality Score Recalibration (BQSR) \n\t\e[94m\e[1m|- \e[0m\e[92mNumber of iterations: ${BOOTS} \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam\e[0m"

# Call SNPs > BQSR > Call SNPs > BQSR > Call SNPs > BQSR > Call SNPs
#cp ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

# Boostrap a number of times
for (( n = 2; n <= ${BOOTS}; n++ ))
do
  echo -e "\e[91m\tProcessing iteration \e[0m\e[96m${n}\e[0m"

  # Sort and Index bootstram bam
  samtools sort ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap
  samtools index ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam



  if [ "${QFIX}" = "T" ]; then

    # GATK UnifiedGenotyper on realn BAM file, qual 20
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_UnifiedGenotyper.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T UnifiedGenotyper \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      --genotyping_mode DISCOVERY \
      --sample_ploidy 1 \
      --min_base_quality_score 20 \
      -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_${n}.vcf \
      -l FATAL \
      --allow_potentially_misencoded_quality_scores \
      2> >(tee "$logfile")

    # Filter SNPs for high quality
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_SNPfilter.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T VariantFiltration \
      -R ${REF} \
      -V ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_${n}.vcf \
      --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
      --filterName "HQ_fail" \
      --genotypeFilterExpression "GQ < 90.0" \
      --genotypeFilterName "GQ_fail" \
      -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_filtered_${n}.vcf \
      -l FATAL \
      --allow_potentially_misencoded_quality_scores \
      2> >(tee "$logfile")

    # Extract SNPs from VCF that passed the filters
    logfile=${OUT}/${SAMPLE}/${SAMPLE}_SelectVariants.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T SelectVariants \
      -R ${REF} \
      --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_filtered_${n}.vcf \
      -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
      -l FATAL \
      --allow_potentially_misencoded_quality_scores \
      --excludeFiltered

    # Analyze patterns of covariation in the sequence dataset
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BaseRecalibrator_pass1.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T BaseRecalibrator \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      -knownSites:VCF ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
      -o ${OUT}/${SAMPLE}/${SAMPLE}_BQSR_${n}.table \
      -l FATAL \
      --allow_potentially_misencoded_quality_scores \
      2> >(tee "$logfile")

    # Do a second pass to analyze covariation post-recalibration
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BaseRecalibrator_pass2.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T BaseRecalibrator \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      -knownSites:VCF ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
      -BQSR ${OUT}/${SAMPLE}/${SAMPLE}_BQSR_${n}.table \
      -o ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_${n}.table \
      -l FATAL \
      --allow_potentially_misencoded_quality_scores \
      2> >(tee "$logfile")

    # Apply the recalibration
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_ApplyRecalibration.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T PrintReads \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      -BQSR ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_${n}.table \
      -o ${OUT}/${SAMPLE}/${SAMPLE}_tmp.bam \
      --emit_original_quals \
     -l FATAL \
      --allow_potentially_misencoded_quality_scores \
      2> >(tee "$logfile")

  fi



  if [ "${QFIX}" = "F" ]; then

    # GATK UnifiedGenotyper on realn BAM file, qual 20
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_UnifiedGenotyper.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T UnifiedGenotyper \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      --genotyping_mode DISCOVERY \
      --sample_ploidy 1 \
      --min_base_quality_score 20 \
      -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_${n}.vcf \
      -l FATAL \
      2> >(tee "$logfile")

    # Filter SNPs for high quality
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_SNPfilter.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T VariantFiltration \
      -R ${REF} \
      -V ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_${n}.vcf \
      --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
      --filterName "HQ_fail" \
      --genotypeFilterExpression "GQ < 90.0" \
      --genotypeFilterName "GQ_fail" \
      -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_filtered_${n}.vcf \
      -l FATAL \
      2> >(tee "$logfile")

    # Extract SNPs from VCF that passed the filters
    logfile=${OUT}/${SAMPLE}/${SAMPLE}_SelectVariants.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T SelectVariants \
      -R ${REF} \
      --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_filtered_${n}.vcf \
      -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
      -l FATAL \
      --excludeFiltered

    # Analyze patterns of covariation in the sequence dataset
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BaseRecalibrator_pass1.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T BaseRecalibrator \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      -knownSites:VCF ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
      -o ${OUT}/${SAMPLE}/${SAMPLE}_BQSR_${n}.table \
      -l FATAL \
      2> >(tee "$logfile")

    # Do a second pass to analyze covariation post-recalibration
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BaseRecalibrator_pass2.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T BaseRecalibrator \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      -knownSites:VCF ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
      -BQSR ${OUT}/${SAMPLE}/${SAMPLE}_BQSR_${n}.table \
      -o ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_${n}.table \
      -l FATAL \
      2> >(tee "$logfile")

    # Apply the recalibration
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_ApplyRecalibration.err
    java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
      -T PrintReads \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      -BQSR ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_${n}.table \
      -o ${OUT}/${SAMPLE}/${SAMPLE}_tmp.bam \
      --emit_original_quals \
     -l FATAL \
      2> >(tee "$logfile")

  fi


  
  cp ${OUT}/${SAMPLE}/${SAMPLE}_tmp.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

done

# Re-index bam file resulting from the above BQSR bootstrapping
samtools index ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

# Remove surplus bootstrap files and realn.bam
rm ${OUT}/${SAMPLE}/${SAMPLE}_tmp.ba*
rm ${OUT}/${SAMPLE}/${SAMPLE}*BQSR*
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn*
rm ${OUT}/${SAMPLE}/${SAMPLE}_realn.b*




#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# Modify environment for Java7
module load bioinfo/Java7 

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# Might need to suppress stdout from programs called within script

echo -e "\n\e[1m\e[34m======================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop \e[0m\e[93mNGS pipeline\e[0m"
echo -e "\e[1m\e[34m======================\e[0m"

# ==============================================================================
# Read in the user passed parameters
# ==============================================================================
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -i <parameters file> -s <sample> [options]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis pipeline has been developed to map Illumina HiSeq reads generated for the SeqApiPop project. All of the parameters in the parameter file require setting. The pipeline depends upon the installation of a range of software: BWA, GATK, Picard, SAMtools, BEDtools, VCFtools and Platypus. In addition, the R packages gsalib, reshape and ggplot2 are required for the before-after BQSR plots. Finally, two Perl scripts are supplied for processing the VCF files output from Platypus and SAMtools.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-i <str>\tPath to input file containing pipeline parameters\e[0m"
  echo -e "\e[92m\t-s <str>\tSample read file name prefix\e[0m"
  echo -e "\e[92m\t-t <int>\tNumber of CPU threads to use subject to software [4]\e[0m"
  echo -e "\e[92m\t-g <int>\tMax amount of RAM (GB) available to Java [48]\e[0m"
  echo -e "\e[92m\t-b <int>\tNumber of bootstrap iterations during BQSR [3]\e[0m"
  echo -e "\n"

}

# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
INFILE=
BOOTS=2
REF=
CHRINFO=
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

while getopts ":i:s:b:" opt; do
  case $opt in
    i) INFILE=${OPTARG};;
    s) SAMPLE=${OPTARG};;
    b) BOOTS=${OPTARG};;
  esac
done

if [[ -z ${INFILE} ]] | [[ -z ${SAMPLE} ]]
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
  if [[ "$name" == "CHRINFO" ]]; then CHRINFO=${value//\"/}; fi
  if [[ "$name" == "IN" ]]; then IN=${value//\"/}; fi
  if [[ "$name" == "OUT" ]]; then OUT=${value//\"/}; fi
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
if [[ "$CHRINFO" == "" ]]; then echo "CHRINFO missing"; exit 0; fi
if [[ "$IN" == "" ]]; then echo "IN path missing"; exit 0; fi
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
echo -e "\e[96mREF: \e[92m${REF}\e[0m"
echo -e "\e[96mCHRINFO: \e[92m${CHRINFO}\e[0m"
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
##echo -e "\e[92m	${OUT}/${SAMPLE}/logs\e[0m"
##mkdir -p ${OUT}/${SAMPLE}/logs
echo -e "\e[92m	${OUT}/${SAMPLE}/metrics\e[0m"
mkdir -p ${OUT}/${SAMPLE}/metrics
echo -e "\e[92m	${OUT}/${SAMPLE}/vcfs\e[0m"
mkdir -p ${OUT}/${SAMPLE}/vcfs


# Number of steps in pipeline, this will need updating if additional steps are added, to provide an idea of how long is remaining whilst running the pipe
STEPS=30



# ==============================================================================
# Mapping
# ==============================================================================

# Map Reads
#READS_1=${SAMPLE}_R1.fastq.gz
#READS_2=${SAMPLE}_R2.fastq.gz
#echo -e "\e[91m[1/${STEPS}] BWA:Map reads [BWA MEM] \n\t\e[94m\e[1m<- \e[0m\e[92m${REF} \n\t\e[94m\e[1m<- \e[0m\e[92m${IN}/${READS_1} \n\t\e[94m\e[1m<- \e[0m\e[92m${IN}/${READS_2} \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}.sam \e[0m"
##${BWA}/bwa mem -M -R @RG"\t"ID:SeqApiPop"\t"SM:${SAMPLE}"\t"PL:illumina"\t"LB:${SAMPLE}"\t"PU:INRA ${REF} ${IN}/${READS_1} ${IN}/${READS_2} > ${OUT}/${SAMPLE}/${SAMPLE}.sam 

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
#java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#  -T RealignerTargetCreator \
#  -R ${REF} \
#  -I ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
#  -o ${OUT}/${SAMPLE}/logs/${SAMPLE}_dup_intervals.list \
#  -l FATAL \
#  2> >(tee "$logfile")

# Perform realignment around INDELs 
#echo -e "\e[91m[6/${STEPS}] GATK:IndelRealigner \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_IndelRealigner.err
#java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#  -T IndelRealigner \
#  -R ${REF} \
#  -I ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
#  -targetIntervals ${OUT}/${SAMPLE}/logs/${SAMPLE}_dup_intervals.list \
#  -o ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
#  -l FATAL \
#  2> >(tee "$logfile")

# Delete surplus files
#echo -e "\e[91m[7/${STEPS}] Removing surplus files\e[0m"
#rm ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam
#rm ${OUT}/${SAMPLE}/${SAMPLE}_dup.bai



# ==============================================================================
# Calculate mapping metrics on realn BAM
# ==============================================================================

# Depth of coverage via GATK
#echo -e "\e[91m[8/${STEPS}] GATK:DepthOfCoverage \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_GATKcov\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_DpethOfCoverage.err
#java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#  -T DepthOfCoverage \
#  -R ${REF} \
#  -I ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
#  -o ${OUT}/${SAMPLE}/metrics/${SAMPLE}_GATKcov \
#  -ct 2 -ct 5 -ct 8 \
#  --omitDepthOutputAtEachBase \
#  --omitIntervalStatistics \
#  --omitLocusTable \
#  -l FATAL \
#  2> >(tee "$logfile")
##${RSCRIPTS}/Plot-GATKcov.R ${SAMPLE} ${OUT}/${SAMPLE}/metrics/${SAMPLE}_GATKcov.sample_summary
##${RSCRIPTS}/Plot-GATKcovStats.R ${SAMPLE} ${OUT}/${SAMPLE}/metrics/${SAMPLE}_GATKcov.sample_statistics

# Flagstat
#echo -e "\e[91m[9/${STEPS}] SAMtools:flagstat \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn.flagstat\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_flagstat.err
#samtools flagstat ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
#  > ${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn.flagstat \
#  2> >(tee "$logfile")
##${RSCRIPTS}/Plot-Flagstat.R ${SAMPLE} ${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn.flagstat

# Alignment Summary Metrics
#echo -e "\e[91m[10/${STEPS}] Picard:CollectAlignmentSummaryMetrics \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_aln.metrics\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardAlignmentMetrics.log
#java -d64 -jar $PICARD/CollectAlignmentSummaryMetrics.jar \
#  INPUT=${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
#  OUTPUT=${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_aln.metrics \
#  REFERENCE_SEQUENCE=${REF} \
#  QUIET=T \
#  VERBOSITY=ERROR \
#  VALIDATION_STRINGENCY=LENIENT \
#  2> >(tee "$logfile")

# Insert Size Metrics
#echo -e "\e[91m[11/${STEPS}] Picard:CollectInsertSizeMetrics \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_ins.metrics \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_ins_hist.pdf\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardInsertSizeMetrics.log
#java -d64 -jar $PICARD/CollectInsertSizeMetrics.jar \
#  INPUT=${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
#  HISTOGRAM_FILE=${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_ins_hist.pdf \
#  METRIC_ACCUMULATION_LEVEL=ALL_READS \
#  OUTPUT=${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_ins.metrics \
#  REFERENCE_SEQUENCE=${REF} \
#  QUIET=T \
#  VERBOSITY=ERROR \
#  VALIDATION_STRINGENCY=LENIENT \
#  2> >(tee "$logfile")

# Quality Score Distribution
#echo -e "\e[91m[12/${STEPS}] Picard:QualityScoreDistribution \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_QSD.metrics \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_QSD.pdf\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardQualityDist.log
#java -d64 -jar $PICARD/QualityScoreDistribution.jar \
#  INPUT=${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
#  CHART_OUTPUT=${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_QSD.pdf \
#  OUTPUT=${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn_QSD.metrics \
#  REFERENCE_SEQUENCE=${REF} \
#  QUIET=T \
#  VERBOSITY=ERROR \
#  VALIDATION_STRINGENCY=LENIENT \
#  2> >(tee "$logfile")

# Breadth of genome coverage
#echo -e "\e[91m[13/${STEPS}] BEDtools:Genomecov \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn.bam.dcov\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_genomecov.err
#bedtools genomecov -ibam ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
#  -g ${CHRINFO} \
#  -max 20 \
#  > ${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn.bam.dcov \
#  2> >(tee "$logfile")

# Poisson distribution relative to genome coverage
#echo -e "\e[91m[14/${STEPS}] R:Coverage-Poisson.R \n\t\e[94m\e[1m<- ${OUT}/${SAMPLE}/metrics/${SAMPLE}_GATKcov.sample_summary \n\t\e[94m\e[1m<- ${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn.bam.dcov \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_Coverage_Poisson.pdf\e[0m"
##${RSCRIPTS}/Coverage-Poisson.R ${OUT} ${SAMPLE} \
#  ${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn.bam.dcov \
#  ${OUT}/${SAMPLE}/metrics/${SAMPLE}_GATKcov.sample_summary
#rm ${OUT}/${SAMPLE}/metrics/${SAMPLE}_realn.bam.dcov.tmp



# ==============================================================================
# BQSR based on boostrapping UnifiedGenotyper and BaseRecalibrator
# ==============================================================================

#echo -e "\e[91m[15/${STEPS}] GATK:Base Quality Score Recalibration (BQSR) \n\t\e[94m\e[1m|- \e[0m\e[92mNumber of iterations: ${BOOTS} \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam\e[0m"

# Call SNPs > BQSR > Call SNPs > BQSR > Call SNPs > BQSR > Call SNPs
#cp ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

# Boostrap a number of times
#for (( n = 1; n <= ${BOOTS}; n++ ))
#do
#  echo -e "\e[91m\tProcessing iteration \e[0m\e[96m${n}\e[0m"

  # Sort and Index bootstram bam
#  samtools sort ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap
#  samtools index ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

  # GATK UnifiedGenotyper on realn BAM file, qual 20
#  logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_UnifiedGenotyper.err
#  java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#    -T UnifiedGenotyper \
#    -R ${REF} \
#    -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
#    --genotyping_mode DISCOVERY \
#    --sample_ploidy 1 \
#    --min_base_quality_score 20 \
#    -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_${n}.vcf \
#    -l FATAL \
#    2> >(tee "$logfile")

  # Filter SNPs for high quality
#  logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_SNPfilter.err
#  java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#    -T VariantFiltration \
#    -R ${REF} \
#    -V ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_${n}.vcf \
#    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
#    --filterName "HQ_fail" \
#    --genotypeFilterExpression "GQ < 90.0" \
#    --genotypeFilterName "GQ_fail" \
#    -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_filtered_${n}.vcf \
#    -l FATAL \
#    2> >(tee "$logfile")

  # Extract SNPs from VCF that passed the filters
#  logfile=${OUT}/${SAMPLE}/${SAMPLE}_SelectVariants.err
#  java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#    -T SelectVariants \
#    -R ${REF} \
#    --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_filtered_${n}.vcf \
#    -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
#    -l FATAL \
#    --excludeFiltered

  # Analyze patterns of covariation in the sequence dataset
#  logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BaseRecalibrator_pass1.err
#  java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#    -T BaseRecalibrator \
#    -R ${REF} \
#    -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
#    -knownSites:VCF ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
#    -o ${OUT}/${SAMPLE}/${SAMPLE}_BQSR_${n}.table \
#    -l FATAL \
#    2> >(tee "$logfile")

  # Do a second pass to analyze covariation post-recalibration
#  logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BaseRecalibrator_pass2.err
#  java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#    -T BaseRecalibrator \
#    -R ${REF} \
#    -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
#    -knownSites:VCF ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
#    -BQSR ${OUT}/${SAMPLE}/${SAMPLE}_BQSR_${n}.table \
#    -o ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_${n}.table \
#    -l FATAL \
#    2> >(tee "$logfile")

  # Apply the recalibration
#  logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_ApplyRecalibration.err
#  java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#    -T PrintReads \
#    -R ${REF} \
#    -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
#    -BQSR ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_${n}.table \
#    -o ${OUT}/${SAMPLE}/${SAMPLE}_tmp.bam \
#    --emit_original_quals \
#    -l FATAL \
#    2> >(tee "$logfile")
  
#  cp ${OUT}/${SAMPLE}/${SAMPLE}_tmp.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

#done

# Cleanup
#rm ${OUT}/${SAMPLE}/${SAMPLE}_tmp.ba*

# Re-index bam file resulting from the above BQSR bootstrapping
#samtools index ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

### Generate before/after plots (R packages gsalib, reshape and ggplot2)
### Currently hashed out because an R update caused problems loading packages
###echo -e "\e[91m[15/${STEPS}] GATK:BQSR before-after plots \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/${SAMPLE}_BQSR_plots_1-${BOOTS}.pdf\e[0m \n\t\e[94m\e[1m** \e[0m\e[93mDisabled due to error in R\e[0m"
###logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_RecalibrationPlots.err
###java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
###  -T AnalyzeCovariates \
###  -R ${REF} \
###  -before ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_1.table \
###  -after ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_${BOOTS}.table \
###  -plots ${OUT}/${SAMPLE}/${SAMPLE}_BQSR_plots_1-${BOOTS}.pdf \
###  2> >(tee "$logfile")

# Quality Score Distribution
#echo -e "\e[91m[16/${STEPS}] Picard:QualityScoreDistribution \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_bootstrap_QSD.metrics \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/metrics/${SAMPLE}_bootstrap_QSD.pdf\e[0m"
#logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardQualityDist_recal.log
#java -d64 -jar $PICARD/QualityScoreDistribution.jar \
#  INPUT=${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
#  CHART_OUTPUT=${OUT}/${SAMPLE}/metrics/${SAMPLE}_bootstrap_QSD.pdf \
#  OUTPUT=${OUT}/${SAMPLE}/metrics/${SAMPLE}_bootstrap_QSD.metrics \
#  REFERENCE_SEQUENCE=${REF} \
#  QUIET=T \
#  VERBOSITY=ERROR \
#  VALIDATION_STRINGENCY=LENIENT \
#  2> >(tee "$logfile")

# Remove surplus bootstrap files and realn.bam
#echo -e "\e[91m[17/${STEPS}] Removing surplus files\e[0m"
#rm ${OUT}/${SAMPLE}/${SAMPLE}*BQSR*
#rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn*
rm ${OUT}/${SAMPLE}/${SAMPLE}_realn.b*



# ==============================================================================
# Perform SNP calling using GATK UnifiedGenotyper, Platypus, SAMtools mpileup
# ==============================================================================

# Platypus
echo -e "\e[91m[18/${STEPS}] Platypus:CallVariants \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_SNPs.vcf\e[0m"
python ${PLATYPUS}/Platypus.py callVariants \
  --output=${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus.vcf \
  --refFile=${REF} \
  --bamFiles=${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
  --nIndividuals=1 \
  --genIndels=1 \
  --genSNPs=1 \
  --minMapQual=30 \
  --minBaseQual=20 \
  --mergeClusteredVariants=1 \
  --logFileName=${OUT}/${SAMPLE}/logs/${SAMPLE}_platypus.log

# Extract homozygous SNPs
echo -e "\e[91m[19/${STEPS}] Perl Script:ParseVCF.pl (Extract homozygous SNPs) \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_SNPs.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_SNPs_hom.vcf\e[0m"
perl ${PERLSCRIPTS}/ParseVCF.pl ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_SNPs_hom.vcf
# Extract SNPs from VCF that passed the filters
logfile=${OUT}/${SAMPLE}/${SAMPLE}_platypus_FinalSelectVariants.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R ${REF} \
  --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_SNPs_hom.vcf \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf \
  -l FATAL \
  --excludeFiltered

# Generate stats
echo -e "\e[91m[20/${STEPS}] VCFtools:vcf-stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.stats\e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.stats
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_SNPs_hom.*
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus*.idx


# GATK UnifiedGenotyper
echo -e "\e[91m[21/${STEPS}] GATK:UnifiedGenotyper \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_SNPs_hom.vcf\e[0m"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_UnifiedGenotyper.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T UnifiedGenotyper \
  -R ${REF} \
  -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
  --genotyping_mode DISCOVERY \
  --genotype_likelihoods_model BOTH \
  --sample_ploidy 1 \
  --output_mode EMIT_VARIANTS_ONLY \
  --min_base_quality_score 20 \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG.vcf \
  -l FATAL \
  2> >(tee "$logfile")

#Filter SNPs for high quality
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_FinalSNPfilter.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R ${REF} \
  -V ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG.vcf \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 30.0" \
  --filterName "HQ_fail" \
  --genotypeFilterExpression "GQ < 50.0" \
  --genotypeFilterName "GQ_fail" \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_filtered.vcf \
  -l FATAL \
  2> >(tee "$logfile")

# Extract SNPs from VCF that passed the filters
logfile=${OUT}/${SAMPLE}/${SAMPLE}_GATK_FinalSelectVariants.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R ${REF} \
  --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_filtered.vcf \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf \
  -l FATAL \
  --excludeFiltered

# Generate stats
echo -e "\e[91m[22/${STEPS}] VCFtools:vcf-stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.stats\e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.stats
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_filtered.* 
rm  ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK*.idx


# SAMtools mpileup on realn BAM file, mapQuality (q) baseQuality (Q)
echo -e "\e[91m[23/${STEPS}] SAMtools:mpileup \n\t\e[94m\e[1m<- \e[0m\e[92m${REF} ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup.vcf\e[0m"
samtools mpileup -C 50 -q 30 -Q 20 -uDf ${REF} ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam | bcftools view -cgv - > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup.vcf

# Extract homozygous Variants
echo -e "\e[91m[24/${STEPS}] Perl Script:ParseVCF.pl (Extract homozygous SNPs) \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_tmp.vcf\e[0m"
perl ${PERLSCRIPTS}/ParseVCF.pl ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_tmp.vcf

# Strip out INDELs
echo -e "\e[91m[25/${STEPS}] Perl Script:PileupStripINDELs.pl (Remove INDELs) \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_tmp.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.vcf\e[0m"
perl ${PERLSCRIPTS}/PileupStripINDELs.pl ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_tmp.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.vcf
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_tmp.vcf

# Generate stats
echo -e "\e[91m[26/${STEPS}] VCFtools:vcf-stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.stats\e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.stats


# Run BAYSIC to select best SNPs
echo -e "\e[91m[27/${STEPS}] BAYSIC:Optimal SNPs \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.vcf \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf\e[0m"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BAYSIC.err
baysic.pl \
  --statsOutFile ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.stats \
  --pvalCutoff 0.8 \
  --vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf \
  --vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.vcf \
  --vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf \
  --countsOutFile ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.counts \
  --vcfOutFile ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf \
  2> >(tee "$logfile")

# Generate stats and Venn Diagram
echo -e "\e[91m[28/${STEPS}] VCFtools:vcf-stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf\e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC_vcf.stats
${RSCRIPTS}/VennDiagram.R ${SAMPLE} ${OUT}/${SAMPLE}/metrics/${SAMPLE}_VCF_VENN.tiff \
  ${SAMPLE} ${OUT}/${SAMPLE}/vcfs/JFM10_TAGCTT_L001_platypus_pass.vcf.gz \
  ${SAMPLE} ${OUT}/${SAMPLE}/vcfs/JFM10_TAGCTT_L001_GATK_UG_pass.vcf.gz \
  ${SAMPLE} ${OUT}/${SAMPLE}/vcfs/JFM10_TAGCTT_L001_pileup_hom.vcf.gz \
  ${SAMPLE} ${OUT}/${SAMPLE}/vcfs/JFM10_TAGCTT_L001_BAYSIC.vcf.gz

# Compress VCF files
bgzip ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus.vcf
bgzip ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf
bgzip ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG.vcf
bgzip ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf
bgzip ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup.vcf
bgzip ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.vcf
bgzip ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf




# ==============================================================================
# SNPs called on MQ 20 according to findings of Ajay et al 2011
# BAYSIC default p-cutoff of 0.8
# ==============================================================================








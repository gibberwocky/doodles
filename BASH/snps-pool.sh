#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# Modify environment for GATK and BCFtools
module load bioinfo/Java7 
module load bioinfo/bcftools

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e


# ==============================================================================
# SNPS.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m======================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop \e[0m\e[93mNGS pipeline\e[0m"
echo -e "\e[1m\e[34m======================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -i <parameters file> -s <sample> [options]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis pipeline has been developed to map Illumina HiSeq reads generated for the SeqApiPop project. All of the parameters in the parameter file require setting. This module requires the installation of GATK, SAMtools, BCFtools, VCFtools, Platypus, BAYSIC, R and the VennDiagrm package.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-i <str>\tPath to input file containing pipeline parameters\e[0m"
  echo -e "\e[92m\t-f <str>\tPath to folder containing fastq.gz files to be mapped\e[0m"
  echo -e "\e[92m\t-s <str>\tSample read file name prefix\e[0m"
  echo -e "\n"

}

# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
INFILE=
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

while getopts ":i:b:n:l:o:" opt; do
  case $opt in
    i) INFILE=${OPTARG};;
    b) BAMS=${OPTARG};;
    n) NAME=${OPTARG};;
    l) LOG=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${INFILE} ]] | [[ -z ${NAME} ]] | [[ -z ${OUT} ]]
then
  usage
  exit 1
fi

# ==============================================================================
# Read in variables
# ==============================================================================

while IFS="=" read name value
do

  if [[ "$name" == "REF" ]]; then REF=${value//\"/}; fi
  if [[ "$name" == "IN" ]]; then IN=${value//\"/}; fi
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

# Create folders for storing files
mkdir -p ${OUT}/${NAME}/logs
mkdir -p ${OUT}/${NAME}/metrics
mkdir -p ${OUT}/${NAME}/vcfs

# Number of steps in pipeline, this will need updating if additional steps are added, to provide an idea of how long is remaining whilst running the pipe
STEPS=11

# ==============================================================================
# Perform SNP calling using GATK UnifiedGenotyper, Platypus, SAMtools mpileup
# ==============================================================================

# Platypus
echo -e "\e[91m[1/${STEPS}] Platypus:CallVariants \e[0m"
python ${PLATYPUS}/Platypus.py callVariants \
  --output=${OUT}/${NAME}_tmp.vcf \
  --refFile=${REF} \
  --bamFiles=${BAMS} \
  --nIndividuals=1 \
  --genIndels=1 \
  --genSNPs=1 \
  --minMapQual=30 \
  --minBaseQual=20 \
  --mergeClusteredVariants=1 \
  --logFileName=${LOG}/${NAME}/${NAME}_platypus.log

# Sort SNPs into correct contig order
DICT=(${REF/.fa/.dict})
perl ${PERLSCRIPTS}/vcfsorter.pl ${DICT} \
  ${OUT}/${NAME}_tmp.vcf > ${OUT}/${NAME}_platypus.vcf

# Extract SNPs from VCF that passed the filters
logfile=${LOG}/${NAME}/${NAME}_platypus_FinalSelectVariants.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R ${REF} \
  --variant ${OUT}_platypus.vcf \
  -o ${OUT}/${NAME}_tmp.vcf \
  -l FATAL \
  --excludeFiltered \
  2> >(tee "$logfile")

# Extract SNPs
bcftools view --types snps -M 2 -O v \
  -o ${OUT}/${NAME}_platypus_pass.vcf  \
  ${OUT}/${NAME}_tmp.vcf

# Generate stats
echo -e "\e[91m[2/${STEPS}] VCFtools:vcf-stats \e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${NAME}_platypus_pass.vcf \
  > ${OUT}/${NAME}_platypus_pass.stats
rm ${OUT}/${NAME}_tmp.vcf


# GATK UnifiedGenotyper
echo -e "\e[91m[3/${STEPS}] GATK:UnifiedGenotyper \e[0m"
logfile=${LOG}/${NAME}/${NAME}_UnifiedGenotyper.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T UnifiedGenotyper \
  -R ${REF} \
  -I ${BAMS} \
  --genotyping_mode DISCOVERY \
  --genotype_likelihoods_model BOTH \
  --output_mode EMIT_VARIANTS_ONLY \
  --min_base_quality_score 20 \
  --max_alternate_alleles 2 \
  -o ${OUT}/${NAME}_GATK_UG.vcf \
  -l FATAL \
  2> >(tee "$logfile")

# Filter SNPs for high quality
logfile=${LOG}/${NAME}/${NAME}_FinalSNPfilter.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R ${REF} \
  -V ${OUT}/${NAME}_GATK_UG.vcf \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 30.0" \
  --filterName "HQ_fail" \
  --genotypeFilterExpression "GQ < 30.0" \
  --genotypeFilterName "GQ_fail" \
  -o ${OUT}/${NAME}_GATK_UG_filtered.vcf \
  -l FATAL \
  2> >(tee "$logfile")

# Extract SNPs from VCF that passed the filters
logfile=${LOG}/${NAME}/${NAME}_GATK_FinalSelectVariants.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R ${REF} \
  --variant ${OUT}/${NAME}_GATK_UG_filtered.vcf \
  -o ${OUT}/${NAME}_tmp.vcf \
  -l FATAL \
  --excludeFiltered

# Extract homozygous SNPs
bcftools view --types snps -M 2 -O v \
  -o ${OUT}/${NAME}_GATK_UG_pass.vcf  \
  ${OUT}/${NAME}_tmp.vcf

# Generate stats
echo -e "\e[91m[4/${STEPS}] VCFtools:vcf-stats \e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${NAME}_GATK_UG_pass.vcf > ${OUT}/${NAME}_GATK_UG_pass.stats
rm ${OUT}/${NAME}_GATK_UG_filtered.* 
rm ${OUT}/${NAME}_GATK*.idx
rm ${OUT}/${NAME}_tmp.*


# SAMtools mpileup on realn BAM file, mapQuality (q) baseQuality (Q)
# Unload BCFtools module to use classic SNP pileup
readarray bam_array < ${BAMS}
module unload bioinfo/bcftools
echo -e "\e[91m[5/${STEPS}] SAMtools:mpileup \e[0m"
samtools mpileup -C 50 -q 30 -Q 20 -uDf ${REF} ${bam_array} | bcftools view -cgv - > ${OUT}/${NAME}_pileup.vcf

# Extract homozygous SNPs
module load bioinfo/bcftools
bcftools view --types snps -M 2 -O v \
  -o ${OUT}/${NAME}_pileup_pass.vcf  \
  ${OUT}/${NAME}_pileup.vcf

# Generate stats
echo -e "\e[91m[6/${STEPS}] VCFtools:vcf-stats \e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${NAME}_pileup_pass.vcf > ${OUT}/${NAME}_pileup_pass.stats

# Run BAYSIC to select best SNPs
echo -e "\e[91m[7/${STEPS}] BAYSIC:Optimal SNPs \e[0m"
logfile=${LOG}/${NAME}/${NAME}_BAYSIC.err
baysic.pl \
  --statsOutFile ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.stats \
  --pvalCutoff 0.8 \
  --vcf ${OUT}/${NAME}_GATK_UG_pass.vcf \
  --vcf ${OUT}/${NAME}_pileup_pass.vcf \
  --vcf ${OUT}/${NAME}_platypus_pass.vcf \
  --countsOutFile ${OUT}/${NAME}_BAYSIC.counts \
  --vcfOutFile ${OUT}/${NAME}_BAYSIC.vcf \
  2> >(tee "$logfile")

# Generate stats and Venn Diagram
echo -e "\e[91m[8/${STEPS}] VCFtools:vcf-stats \e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${NAME}_BAYSIC.vcf > ${OUT}/${NAME}_BAYSIC_vcf.stats
${RSCRIPTS}/VennDiagram.R ${NAME} ${OUT}/${NAME}_VCF_VENN.tiff \
  ${OUT}/${NAME}_platypus_pass.vcf \
  ${OUT}/${NAME}_GATK_UG_pass.vcf \
  ${OUT}/${NAME}_pileup_pass.vcf \
  ${OUT}/${NAME}_BAYSIC.vcf

# Compress VCF files
bgzip -f ${OUT}/${NAME}_platypus.vcf
bgzip -f ${OUT}/${NAME}_platypus_pass.vcf
bgzip -f ${OUT}/${NAME}_GATK_UG.vcf
bgzip -f ${OUT}/${NAME}_GATK_UG_pass.vcf
bgzip -f ${OUT}/${NAME}_pileup.vcf
bgzip -f ${OUT}/${NAME}_pileup_pass.vcf
bgzip -f ${OUT}/${NAME}_BAYSIC.vcf




#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# Modify environment
module load bioinfo/Java7 
module load bioinfo/bcftools

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
  echo -e "\e[92m\t-f <str>\tPath to folder containing fastq.gz files to be mapped\e[0m"
  echo -e "\e[92m\t-s <str>\tSample read file name prefix\e[0m"
  echo -e "\e[92m\t-b <int>\tNumber of bootstrap iterations during BQSR [3]\e[0m"
  echo -e "\n"

}

# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
INFILE=
SITES=
OUTFILE=
BAMS=

while getopts ":i:b:s:o:" opt; do
  case $opt in
    i) INFILE=${OPTARG};;
    b) BAMS=${OPTARG};;
    s) SITES=${OPTARG};;
    o) OUTFILE=${OPTARG};;
  esac
done

if [[ -z ${INFILE} ]] | [[ -z ${SITES} ]] | [[ -z ${OUTFILE} ]] | [[ -z ${BAMS} ]]
then
  exit 1
fi


# ==============================================================================
# Read in variables
# ==============================================================================

while IFS="=" read name value
do

  if [[ "$name" == "REF" ]]; then REF=${value//\"/}; fi
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
if [[ "$BWA" == "" ]]; then echo "BWA missing"; exit 0; fi
if [[ "$PICARD" == "" ]]; then echo "PICARD path missing"; exit 0; fi
if [[ "$GATK" == "" ]]; then echo "GATK path missing"; exit 0; fi
if [[ "$PLATYPUS" == "" ]]; then echo "PLATYPUS path missing"; exit 0; fi
if [[ "$VCFT" == "" ]]; then echo "VCFT path missing"; exit 0; fi
if [[ "$BAYSIC" == "" ]]; then echo "BAYSIC path missing"; exit 0; fi
if [[ "$PERLSCRIPTS" == "" ]]; then echo "PERLSCRIPTS path missing"; exit 0; fi
if [[ "$RSCRIPTS" == "" ]]; then echo "RSCRIPTS path missing"; exit 0; fi

echo -e "\n\e[91mVariables imported succesfully:\e[0m"

echo -e "\e[96mREF: \e[92m${REF}\e[0m"
echo -e "\e[96mBWA: \e[92m${BWA}\e[0m"
echo -e "\e[96mPICARD: \e[92m${PICARD}\e[0m"
echo -e "\e[96mGATK: \e[92m${GATK}\e[0m"
echo -e "\e[96mPLATYPUS: \e[92m${PLATYPUS}\e[0m"
echo -e "\e[96mVCFT: \e[92m${VCFT}\e[0m"
echo -e "\e[96mBAYSIC: \e[92m${BAYSIC}\e[0m"
echo -e "\e[96mPERLSCRIPTS: \e[92m${PERLSCRIPTS}\e[0m"
echo -e "\e[96mRSCRIPTS: \e[92m${RSCRIPTS}\e[0m"


# ====================================================================
# Prior to running this script:
# ====================================================================

# 1) Merge sample VCFs from a single population into 1 file
#	POP="JFM"
#	VCFS=(/home/dwragg/work/tmp/*/vcfs/${POP}*clean.vcf.gz)
#	VCFOUT="${DUMP}/vcfs/${POP}_sites.vcf.gz"
# 2) Generate VCF file containing merged sites but not genotypes
#	bcftools merge ${VCFS[@]} | bcftools view --types snps -G -M 2 -g ^het -O z -o ${VCFOUT} -
#	tabix -p vcf ${VCFOUT}
# 3) Generate list of BAM files of individuals
#	BAMS=(/home/dwragg/work/tmp/${POP}*/${POP}*.bam)
#	BAMOUT="${DUMP}/vcfs/${POP}_BAM.list"
#	ls ${BAMS[@]} > ${BAMOUT}

# ====================================================================

# Genotype samples with GATK using sites from previous step
logfile=${OUTFILE}.err
java -d64 -Xmx48g -jar ${GATK}/GenomeAnalysisTK.jar \
  -T UnifiedGenotyper \
  -R ${REF} \
  -I ${BAMS} \
  --alleles ${SITES} \
  -L ${SITES} \
  --genotyping_mode GENOTYPE_GIVEN_ALLELES \
  --genotype_likelihoods_model SNP \
  --sample_ploidy 2 \
  --output_mode EMIT_ALL_SITES \
  --min_base_quality_score 20 \
  --standard_min_confidence_threshold_for_calling 50.0 \
  -G none \
  -o ${OUTFILE} \
  2> >(tee "$logfile")



#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# Load environment
module load bioinfo/bcftools
module load bioinfo/Java7



# ==============================================================================
# IMPUTE.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m======================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop \e[0m\e[93mNGS pipeline\e[0m"
echo -e "\e[1m\e[34m======================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -i <samples file> -o <output path> -p <script path> -n <output name>\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis pipeline has been developed to map Illumina HiSeq reads generated for the SeqApiPop project. This module merges the VCF files of samples from different populations, and uses a set of master sites for genotyping each population independently. Missing genotypes are imputed for each population independently using Beagle. After, the resulting population VCFs are merged into a single VCF file.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-i <str>\tPath to tab-delimited input file containing population assignment [int], path to VCF file, path to BAM file, for each sample\e[0m"
  echo -e "\e[92m\t-o <str>\tPath to folder in which to write the resulting files\e[0m"
  echo -e "\e[92m\t-p <str>\tPath to folder containing Perl script SubHet4Miss.pl\e[0m"
  echo -e "\e[92m\t-n <int>\tNaming prefix for final merged vCF file\e[0m"
  echo -e "\n"

}




# ==============================================================================
# IMPUTE
# ==============================================================================

BEAGLE=/usr/local/bioinfo/src/Beagle/beagle4.jar
PERLSCRIPTS=/home/dwragg/save/Scripts/Perl
PIPE=
GATK=/usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.1-1
REF=/home/dwragg/save/Apis/Apis_mellifera.fa
INFILE=
DUMP=
BAMS=
VCFS=
HAPF=
POPS=
MAX=
MERGER=
NAME=
SITES=

while getopts ":i:o:p:g:r:x:z:n:s:" opt; do
  case $opt in
    i) INFILE=${OPTARG};;
    o) DUMP=${OPTARG};;
    p) PIPE=${OPTARG};;
    g) GATK=${OPTARG};;
    r) REF=${OPTARG};;
    x) PERLSCRIPTS=${OPTARG};;
    z) BEAGLE=${OPTARG};;
    n) NAME=${OPTARG};;
    s) SITES=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${DUMP} ]] | [[ -z ${PIPE} ]] | [[ -z ${REF} ]] | [[ -z ${GATK} ]] | [[ -x ${PERLSCRIPTS} ]] | [[ -z ${BEAGLE} ]] | [[ -z ${NAME} ]] | [[ -z ${SITES} ]]
then
  usage
  exit 1
fi

# ==============================================================================
# Read in variables
# ==============================================================================

POPS=($(cut -f1 ${INFILE})) # list of pops
VCFS=($(cut -f2 ${INFILE})) # list of vcfs
BAMS=($(cut -f3 ${INFILE})) # list of bams
HAPF=($(cut -f4 ${INFILE})) # list of bams
IFS=$'\n'; POPS_MAX=($(echo "${POPS[*]}" | sort -nr | head -n1)) # max pops



# Iterate through each population
POP=1
for ((POP=1; POP <= ${POPS_MAX}; ++POP)); do
  POS_START="-1"
  for ((i=0; i < ${#POPS[@]}; ++i)); do
	if [ "${POP}" = "${POPS[$i]}" ]; then
		if [ "${POS_START}" = "-1" ]; then 
			POS_START=${i}
		fi
        	POS_END=$((${i}+1))
	fi
  done
  POS_END=$((${POS_END}-${POS_START}))

  # Generate file containing list of BAM files to operate on
  BAMFILE=${DUMP}/${POP}_BAMS.list
  printf "%s\n" "${BAMS[@]:${POS_START}:${POS_END}}" > ${BAMFILE}

printf "\nPOP:\t%s" "${POP}"
#printf "\nFROM:\t%s" "${POS_START}"
#printf "\nTO:\t%s" "${POS_END}"
#POS_END=$((${POS_END}-${POS_START}))
printf "\nN:\t%s\n" ${POS_END}
printf "%s\n" "${BAMS[@]:${POS_START}:${POS_END}}"


# if [[ ${POP} > 1 ]]; then

  # Genotype samples with GATK using sites from previous step (long-time)
  echo -e "\e[91mGenotyping samples from population \e[94m${POP} \e[91mwith GATK UnifiedGenotyper \e[0m"
  cd ${DUMP}
  java -d64 -Xmx48g -jar ${GATK}/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R ${REF} \
    -I ${BAMFILE} \
    --alleles ${SITES} \
    -L ${SITES} \
    --genotyping_mode GENOTYPE_GIVEN_ALLELES \
    --genotype_likelihoods_model SNP \
    --sample_ploidy 2 \
    --output_mode EMIT_ALL_SITES \
    --min_base_quality_score 20 \
    --standard_min_confidence_threshold_for_calling 50.0 \
    -G none \
    -o ${DUMP}/${POP}_total.vcf.gz

# fi



  # Perl script to replace heterozygous genotypes as missing if haploid
  if [ "${HAPF[${POP_START}]}" = "Y" ]; then
    echo -e "\e[91mExecute SubHet4Miss.pl to replace heterozygous genotypes as missing in population \e[94m${POP} \e[0m"
    gunzip ${DUMP}/${POP}_total.vcf.gz
    perl ${PERLSCRIPTS}/SubHet4Miss.pl ${DUMP}/${POP}_total.vcf >  ${DUMP}/${POP}_subhet.vcf
    rm ${DUMP}/${POP}_total.*
    bgzip ${DUMP}/${POP}_subhet.vcf
    tabix -p vcf ${DUMP}/${POP}_subhet.vcf.gz
  fi
  # Otherwise just rename file
  if [ "${HAPF[${POP_START}]}" = "N" ]; then
    mv ${DUMP}/${POP}_total.vcf.gz ${DUMP}/${POP}_subhet.vcf.gz
    tabix -p vcf ${DUMP}/${POP}_subhet.vcf.gz
  fi

  # Filter SNPs to remove SNPs with no calls (~2.5% per populations)
  echo -e "\e[91mFilter SNPs to remove those with 0% call rate, typically 2.5%, in population \e[94m${POP} \e[0m"
  bcftools view --types snps -U -M 2 -O z \
    -o ${DUMP}/${POP}_called.vcf.gz ${DUMP}/${POP}_subhet.vcf.gz
  tabix -fp vcf ${DUMP}/${POP}_called.vcf.gz

  # Write uncalled sites to file (may want them at some point?)
  echo -e "\e[91mRecord sites with 0% call rate to file for population \e[94m${POP} \e[0m"
  bcftools view --types snps -G -u -O z \
    -o ${DUMP}/${POP}_uncalled.vcf.gz ${DUMP}/${POP}_subhet.vcf.gz
  rm ${DUMP}/${POP}_subhet.*

  # Impute missing data if haploid
  if [ "${HAPF[${POP_START}]}" = "Y" ]; then
    echo -e "\e[91mImpute missing genotypes with BEAGLE in population \e[94m${POP} \e[0m"
    cd ${DUMP}
    java -d64 -jar ${BEAGLE} \
      gt=${DUMP}/${POP}_called.vcf.gz \
      out=${DUMP}/${POP}_beagle \
      usephase=true
    tabix -fp vcf ${DUMP}/${POP}_beagle.vcf.gz
    rm ${DUMP}/${POP}_called.*
  fi
  # Phase and impute missing data if diploid
  if [ "${HAPF[${POP_START}]}" = "N" ]; then
    echo -e "\e[91mImpute missing genotypes with BEAGLE in population \e[94m${POP} \e[0m"
    cd ${DUMP}
    java -d64 -jar ${BEAGLE} \
      gt=${DUMP}/${POP}_called.vcf.gz \
      out=${DUMP}/${POP}_beagle \
      usephase=false
    tabix -fp vcf ${DUMP}/${POP}_beagle.vcf.gz
    rm ${DUMP}/${POP}_called.*
  fi

  # Record for merger
  MERGER=("${MERGER[@]}" "${DUMP}/${POP}_beagle.vcf.gz")

done

# Merge VCFs from each population
echo -e "\e[91mMerge population data into a single VCF, \e[94m${NAME}.vcf.gz \e[0m"
printf "\t\e[94m<- \e[0m\e[96m%s\n\e[0m" "${MERGER[@]:1:${#MERGER[@]}}"
bcftools merge ${MERGER[@]} > ${DUMP}/${NAME}.vcf
bgzip ${DUMP}/${NAME}.vcf
tabix -fp vcf ${DUMP}/${NAME}.vcf.gz
rm ${MERGER[@]}






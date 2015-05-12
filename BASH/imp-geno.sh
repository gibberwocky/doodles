#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# Load environment
module load bioinfo/bcftools
module load bioinfo/Java7


# ==============================================================================
# imp-geno.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Genotype / impute sites\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -o <output path> -s <VCF file of sites> -v <file listing VCF files> -b <file listing BAM files> -n <output name>\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mGenotypes a population of samples for a given set of sites and imputes missing genotypes within the population.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-i <str>\tPath to input file containing pipeline parameters [defafult = /save/seqapipop/Scripts/params]\e[0m"
  echo -e "\e[92m\t-p <int>\tPloidy [default = 1] determines whether or not to substitute heterozygous genotypes for missing\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mGATK\e[0m"
  echo -e "\e[92mSubHet4Miss.pl\e[0m"
  echo -e "\e[92mbcftools\e[0m"
  echo -e "\e[92mBeagle\e[0m"
  echo -e "\n"

}




# ==============================================================================
# IMPUTE
# ==============================================================================

INFILE=/save/seqapipop/Scripts/params
PIPE=/save/seqapipop/Scripts
VCFS=
BAMS=
NAME=
SITES=
PLOID=1
DUMP=
BWA=
PICARD=
GATK=
PLATYPUS=
VCFT=
PERLSCRIPTS=
RSCRIPTS=
BAYSIC=
BEAGLE=

while getopts ":i:o:v:b:n:s:p:" opt; do
  case $opt in
    i) INFILE=${OPTARG};;
    o) DUMP=${OPTARG};;
    v) VCFS=${OPTARG};;
    b) BAMS=${OPTARG};;
    n) NAME=${OPTARG};;
    s) SITES=${OPTARG};;
    p) PLOID=${OPTARG};;
  esac
done

if [[ -z ${INFILE} ]] | [[ -z ${DUMP} ]] | [[ -z ${VCFS} ]] | [[ -z ${BAMS} ]] | [[ -z ${NAME} ]] | [[ -z ${SITES} ]]
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


# Genotype samples with GATK using sites from previous step (long-time)
echo -e "\e[91mGenotyping samples from population \e[94m${NAME} \e[91mwith GATK UnifiedGenotyper \e[0m"
cd ${DUMP}
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
  --standard_min_confidence_threshold_for_calling 30.0 \
  -G none \
  -o ${DUMP}/${NAME}.vcf.gz



# Perl script to replace heterozygous genotypes as missing if haploid
 if [ "${PLOID}" = "1" ]; then
    echo -e "\e[91mExecute SubHet4Miss.pl to replace heterozygous genotypes as missing in population \e[94m${NAME} \e[0m"
    gunzip ${DUMP}/${NAME}.vcf.gz
    perl ${PERLSCRIPTS}/SubHet4Miss.pl ${DUMP}/${NAME}.vcf \
      > ${DUMP}/${NAME}_subhet.vcf
    rm ${DUMP}/${NAME}.vcf
    bgzip ${DUMP}/${NAME}_subhet.vcf
 fi
# Otherwise just rename file
 if [ "${PLOID}" = "2" ]; then
    mv ${DUMP}/${NAME}.vcf.gz ${DUMP}/${NAME}_subhet.vcf.gz
fi
tabix -p vcf ${DUMP}/${NAME}_subhet.vcf.gz

# Filter SNPs to remove SNPs with no calls (~2.5% per populations)
echo -e "\e[91mFilter SNPs to remove those with 0% call rate, typically 2.5%, in population \e[94m${NAME} \e[0m"
bcftools view --types snps -U -M 2 -O z \
  -o ${DUMP}/${NAME}_preImp.vcf.gz ${DUMP}/${NAME}_subhet.vcf.gz
tabix -fp vcf ${DUMP}/${NAME}_preImp.vcf.gz

# Impute missing data
echo -e "\e[91mImpute missing genotypes with BEAGLE in population \e[94m${NAME} \e[0m"
cd ${DUMP}
java -d64 -jar ${BEAGLE} \
  gt=${DUMP}/${NAME}_preImp.vcf.gz \
  out=${DUMP}/${NAME}_beagle
tabix -fp vcf ${DUMP}/${NAME}_beagle.vcf.gz
#rm ${DUMP}/${NAME}_preImp.*

# Clean up
#rm ${DUMP}/${NAME}_subhet*
#rm ${DUMP}/${NAME}.vcf.gz.idx




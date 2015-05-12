#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# Modify environment for GATK and BCFtools
module load bioinfo/Java7 
module load bioinfo/bcftools

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e


# ==============================================================================
# snps-diploid.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Diploid SNP calling\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -s <sample> -f <path to fastq reads> -o <path to output> [options]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis pipeline has been developed to calls SNPs in diploid samples. All of the parameters in the parameter file required to be set"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-i <str>\tPath to input file containing pipeline parameters [defafult = /save/seqapipop/Scripts/params]\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mGATK\e[0m"
  echo -e "\e[92mSAMtools\e[0m"
  echo -e "\e[92mBCFtools\e[0m"
  echo -e "\e[92mVCFtools\e[0m"
  echo -e "\e[92mPlatypus\e[0m"
  echo -e "\e[92mBAYSIC\e[0m"
  echo -e "\e[92mR VennDiagram package\e[0m"
  echo -e "\e[92mPvcfsorter.pl\e[0m"
  echo -e "\n"
}

# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
INFILE=/save/seqapipop/Scripts/params
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
BEAGLE=

while getopts ":i:f:s:o:" opt; do
  case $opt in
    i) INFILE=${OPTARG};;
    f) IN=${OPTARG};;
    s) SAMPLE=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${INFILE} ]] | [[ -z ${SAMPLE} ]] | [[ -z ${OUT} ]]
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
  if [[ "$name" == "BEAGLE" ]]; then BEAGLE=${value//\"/}; fi

done < ${INFILE}

if [[ "$REF" == "" ]]; then echo "REF missing"; exit 0; fi
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

# Create folders for storing files
mkdir -p ${OUT}/${SAMPLE}/logs
mkdir -p ${OUT}/${SAMPLE}/metrics
mkdir -p ${OUT}/${SAMPLE}/vcfs

# Number of steps in pipeline, this will need updating if additional steps are added, to provide an idea of how long is remaining whilst running the pipe
STEPS=11

# ==============================================================================
# Perform SNP calling using GATK UnifiedGenotyper, Platypus, SAMtools mpileup
# ==============================================================================

# Platypus
echo -e "\e[91m[1/${STEPS}] Platypus:CallVariants \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus.vcf\e[0m"
python ${PLATYPUS}/Platypus.py callVariants \
  --output=${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf \
  --refFile=${REF} \
  --bamFiles=${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
  --nIndividuals=1 \
  --genIndels=1 \
  --genSNPs=1 \
  --minMapQual=30 \
  --minBaseQual=20 \
  --mergeClusteredVariants=1 \
  --logFileName=${OUT}/${SAMPLE}/logs/${SAMPLE}_platypus.log

# Sort SNPs into correct contig order
DICT=(${REF/.fa/.dict})
perl ${PERLSCRIPTS}/vcfsorter.pl ${DICT} \
  ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus.vcf

# Extract SNPs from VCF that passed the filters
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_platypus_FinalSelectVariants.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T SelectVariants \
  -R ${REF} \
  --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus.vcf \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf \
  -l FATAL \
  --excludeFiltered \
  2> >(tee "$logfile")

# Extract SNPs with 2 alleles
bcftools view --types snps -M 2 -O v \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf  \
  ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf

# Generate stats
echo -e "\e[91m[2/${STEPS}] VCFtools:vcf-stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.stats\e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.stats
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.*


# GATK UnifiedGenotyper
echo -e "\e[91m[3/${STEPS}] GATK:UnifiedGenotyper \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_SNPs_hom.vcf\e[0m"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_UnifiedGenotyper.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T UnifiedGenotyper \
  -R ${REF} \
  -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
  --genotyping_mode DISCOVERY \
  --genotype_likelihoods_model BOTH \
  --output_mode EMIT_VARIANTS_ONLY \
  --min_base_quality_score 20 \
  --max_alternate_alleles 2 \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG.vcf \
  -l FATAL \
  2> >(tee "$logfile")

# Filter SNPs for high quality
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_FinalSNPfilter.err
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R ${REF} \
  -V ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG.vcf \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 30.0" \
  --filterName "HQ_fail" \
  --genotypeFilterExpression "GQ < 30.0" \
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
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf \
  -l FATAL \
  --excludeFiltered

# Extract SNPs with 2 alleles
bcftools view --types snps -M 2 -O v \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf  \
  ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf

# Generate stats
echo -e "\e[91m[4/${STEPS}] VCFtools:vcf-stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.stats\e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.stats
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_filtered.* 
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK*.idx
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.*


# SAMtools mpileup on realn BAM file, mapQuality (q) baseQuality (Q)
echo -e "\e[91m[5/${STEPS}] SAMtools:mpileup \n\t\e[94m\e[1m<- \e[0m\e[92m${REF} ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup.vcf\e[0m"
samtools mpileup -uBg -f ${REF} -C 50 -q 30 -Q 20 -t DP ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam | bcftools call -vmO v -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup.vcf

# Extract SNPs with 2 alleles
bcftools view --types snps -M 2 -O v \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf  \
  ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup.vcf

# Generate stats
echo -e "\e[91m[6/${STEPS}] VCFtools:vcf-stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.stats\e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.stats

# Run BAYSIC to select best SNPs
echo -e "\e[91m[7/${STEPS}] BAYSIC:Optimal SNPs \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_hom.vcf \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf\e[0m"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BAYSIC.err
baysic.pl \
  --statsOutFile ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.stats \
  --pvalCutoff 0.8 \
  --vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf \
  --vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf \
  --vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf \
  --countsOutFile ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.counts \
  --vcfOutFile ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf \
  2> >(tee "$logfile")

# Generate stats and Venn Diagram
echo -e "\e[91m[8/${STEPS}] VCFtools:vcf-stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf\e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC_vcf.stats
${RSCRIPTS}/VennDiagram.R ${SAMPLE} ${OUT}/${SAMPLE}/metrics/${SAMPLE}_VCF_VENN.tiff \
  ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf \
  ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf \
  ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf \
  ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf


# Fix 4.2-isms in mpileup for GATK Combine Variants
sed -i 's/VCFv4.2/VCFv4.1/' ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf
sed -i 's/,Version=3>/>/' ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf
sed -i 's/,Version=\"3\">/>/' ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf
sed -i 's/Number=R/Number=./' ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf


# Compress VCF files
bgzip -f ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus.vcf
bgzip -f ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf
bgzip -f ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG.vcf
bgzip -f ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf
bgzip -f ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup.vcf
bgzip -f ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf
bgzip -f ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf



# ==============================================================================
# Generate a single column VCF for the sample
# ==============================================================================


# Index VCFs
tabix -fp vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf.gz
tabix -fp vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf.gz 
tabix -fp vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf.gz 

# Prep interval list from BAYSIC for GATK
cp ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.vcf.pos ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.intervals
sed -i '1d' ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.intervals
sed -i 's/\t/:/g' ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.intervals

# CombineVariants using BAYSIC intervals list
echo -e "\e[91m[9/${STEPS}] GATK:CombineVariants \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf.gz \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_pileup_pass.vcf.gz \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_platypus.vcf.gz \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf\e[0m"
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T CombineVariants \
  -R ${REF} \
  --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_UG_pass.vcf.gz \
  --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_pileup_pass.vcf.gz \
  --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_platypus_pass.vcf.gz \
  -L ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_BAYSIC.intervals \
  -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf \
  --assumeIdenticalSamples \
  --genotypemergeoption UNSORTED \
  -l FATAL

# Use vcf-merge to strip out duplicate sites
echo -e "\e[91m[10/${STEPS}] VCTtools:vcf-merge -d -s \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf.gz \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vcf\e[0m"
bgzip -f ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf
tabix -fp vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf.gz
perl ${VCFT}/vcf-merge -d -s ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.vcf.gz > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vcf
bgzip -f ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vcf
tabix -fp vcf ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vcf.gz

# Generate stats
echo -e "\e[91m[11/${STEPS}] VCTtools:vcf-stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vcf.gz \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.stats\e[0m"
perl ${VCFT}/vcf-stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vcf.gz > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.stats

# Calculate ts/tv rate
echo -e "\e[91m[11/${STEPS}] BCFtools:stats \n\t\e[94m\e[1m<- \e[0m\e[92m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vcf.gz \n\t\e[94m\e[1m-> \e[0m\e[96m${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vchk\e[0m"
bcftools stats ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vcf.gz > ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_clean.vchk


# Clean up
set +e
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_tmp.*
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}*platypus*gz
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}*pileup*gz
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}*GATK*gz
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}*BAYSIC*gz
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}*platypus*gz.tbi
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}*pileup*gz.tbi
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}*GATK*gz.tbi




# ==============================================================================
# SNPs called on MQ 20 according to findings of Ajay et al 2011
# BAYSIC default p-cutoff of 0.8
# ==============================================================================









#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# Modify environment for Java7 and bcftools
module load bioinfo/Java7 
module load bioinfo/bcftools


# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# GC-content.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: GC content by interval\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"

# Clear variables 
PIPE=/save/seqapipop/Scripts
INFILE=/save/seqapipop/Scripts/params
REF=/save/dwragg/Apis/ensembl/Apis_mellifera.GCA_000002195.1.25.dna.genome.fa
DICT=/save/dwragg/Apis/ensembl/Apis_mellifera.GCA_000002195.1.25.dna.genome.dict
OUT=
SAMPLE=
WIN=5000
INC=1000
VCF=

while getopts ":w:i:o:s:v:" opt; do
 case $opt in
   i) INC=${OPTARG};;
   w) WIN=${OPTARG};;
   o) OUT=${OPTARG};;
   s) SAMPLE=${OPTARG};;
   v) VCF=${OPTARG};;
 esac
done

if [[ -z ${OUT} ]] | [[ -z ${INTERVALS} ]] | [[ -z ${VCF} ]] | [[ -z ${SAMPLE} ]]
then
 exit 1
fi

# Read required variables from params file
while IFS="=" read name value
do
 if [[ "$name" == "GATK" ]]; then GATK=${value//\"/}; fi
 if [[ "$name" == "PICARD" ]]; then PICARD=${value//\"/}; fi
done < ${INFILE}



# Extract sample from VCF
bcftools view --types snps -M 2 -O v -s ${SAMPLE} -o ${OUT}/${SAMPLE}_tmp.vcf ${VCF}

# Re-order contigs to match reference
perl ${PIPE}/vcfsorter.pl ${DICT} ${OUT}/${SAMPLE}_tmp.vcf > ${OUT}/${SAMPLE}.vcf
rm ${OUT}/${SAMPLE}_tmp.vcf

# Generate fasta
samtools faidx ${REF} '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' | bcftools consensus ${VCF} -s ${SAMPLE} -o ${OUT}/${SAMPLE}.fa

# Calculate GC in windows
perl ${PIPE}/Perl/calculateGC.pl --fasta ${OUT}/${SAMPLE}.fa --winlen ${WIN} --increment ${INC} > ${OUT}/${SAMPLE}.gc

# Clean-up
rm ${OUT}/${SAMPLE}.vcf
rm ${OUT}/${SAMPLE}.fa





# OLD CODE THAT IS MUCH SLOWER



# Switch to new reference
#REF=${OUT}/${SAMPLE}

# Index reference and create sequence dictionary
#samtools faidx ${REF}.fa
#java -d64 -Xmx48g -jar ${PICARD}/CreateSequenceDictionary.jar \
#  REFERENCE=${REF}.fa \
#  OUTPUT=${REF}.dict

# GC content walker
#ints=($(cat ${INTERVALS}))
#i=0
#for interval in ${ints[@]}
#do
#  i=$((i+1))
#  java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
#   -T GCContentByInterval \
#   -R ${REF}.fa \
#   -o ${OUT}/${SAMPLE}_${i}.tmp \
#   -L ${interval} \
#   -l FATAL
#done
#cat ${OUT}/${SAMPLE}_*.tmp > ${OUT}/${SAMPLE}_GC
#rm ${OUT}/${SAMPLE}*.tmp


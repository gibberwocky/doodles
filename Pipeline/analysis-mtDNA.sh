#!/bin/bash
# IEKNGHER1ujm

# Load modules
module load bioinfo/bcftools

# Paths to main folders
DUMP=/home/dwragg/work/Analysis
PIPE=/save/seqapipop/Scripts
PERL=${PIPE}/Perl
SEQA=/save/seqapipop/Data/Apis-mellifera
AMEL=/work/dwragg/Apis-mellifera

# Merged master sites: 		HARP-JFM-OTH.vcf.gz
# Harpur_preImp.vcf.gz wc -l = 10950024 (header = 5338)
# JFMOTH_preImp.vcf.gz wc -l = 10981674 (header = 5338)

# This is combined Harpur and JFM-OTH SNPs post-filtered for depth of coverage

# Output path and master site file
OUT=${DUMP}/vcfs/HARP-JFM-OTH-15012015
sites=${DUMP}/sites/HARP-JFM-OTH.vcf.gz

# As not imputing, can pool samples into 2 runs: (1) haploid and (2) diploid
Harpur=(${VCF_A[@]} ${VCF_M[@]} ${VCF_C[@]} ${VCF_O[@]})
printf "%s\n" "${Harpur[@]}" > ${DUMP}/vcfs/lists/Harpur_vcfs.list
Harpur=(${BAM_A[@]} ${BAM_M[@]} ${BAM_C[@]} ${BAM_O[@]})
printf "%s\n" "${Harpur[@]}" > ${DUMP}/vcfs/lists/Harpur_bams.list


# 1) extract mtDNA consensus sequence
# 2) generate consensus alignment
REF=/home/dwragg/save/Apis/Apis_mellifera.fa
REGION="NC_001566.1"

A=${DUMP}/vcfs/lists/A_Harpur_samples.list
M=${DUMP}/vcfs/lists/M_Harpur_samples.list
C=${DUMP}/vcfs/lists/C_Harpur_samples.list
O=${DUMP}/vcfs/lists/O_Harpur_samples.list
IDS=( $( cat ${A} ${M} ${C} ${O} ) )
VCFS=/home/dwragg/work/Apis-mellifera
for ID in ${IDS[@]}
do
  samtools faidx ${REF} ${REGION} | \
    bcftools consensus ${VCFS}/*/*/${ID}_clean.vcf.gz \
    -s ${ID} -o ${DUMP}/mtDNA/${ID}_mt.fa
  sed -i "s/^>.*/>${ID}/1" ${DUMP}/mtDNA/${ID}_mt.fa
done

JFMOTH=${DUMP}/vcfs/lists/JFMOTH_samples.list
IDS=($( cat ${JFMOTH} ) )
VCFS=/work/dwragg/Apis-mellifera/SeqApiPop
for ID in ${IDS[@]}
do
  samtools faidx ${REF} ${REGION} | \
    bcftools consensus ${VCFS}/*/*/${ID}_clean.vcf.gz \
    -s ${ID} -o ${DUMP}/mtDNA/${ID}_mt.fa
  sed -i "s/^>.*/>${ID}/1" ${DUMP}/mtDNA/${ID}_mt.fa
done


# Concatenate individual sequences into single file
cat ${DUMP}/mtDNA/*_mt.fa > ${DUMP}/mtDNA/all_mt.FA


# Align sequences and run JModelTest2, check script for -quicktree
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
  -o ${DUMP}/mtDNA \
  -e ${DUMP}/mtDNA \
  ${PIPE}/phylo.sh \
    -i ${DUMP}/mtDNA/all_mt.FA \
    -o ${DUMP}/mtDNA/all_mt.nxs

# Before running MrBayes, check output of jModeltest *.nxs.jmt
# statefreqr correspond to jModelTest f(a) f(c) f(g) f(t)
# revmatpr correspond to jModelTest Ra Rb Rc Rd Re Rf
# MrBayes
module load bioinfo/MrBayes_3.2.2
MrBayes=/usr/local/bioinfo/src/MrBayes/current
cd ${DUMP}/mtDNA
${MrBayes}/mb 
execute all_mt.nxs

# Parameter set 1
lset nucmodel=4by4 nst=6 rates=equal
prset statefreqpr = fixed(0.43, 0.1, 0.05, 0.42)
prset revmatpr=fixed(2.140, 27.536, 2.300, 0.940, 46.066,1.000)

# Parameter set 2
lset nucmodel=4by4 nst=6 rates=gamma

# Run and calculate tree
mcmc ngen=1000000 samplefreq=100 printfreq=10000 diagnfreq=10000 checkfreq=100000
sump
sumt

















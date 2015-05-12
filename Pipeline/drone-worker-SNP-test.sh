#!bin/bash
# IEKNGHER1ujm

# Modify environment for Java7
module load bioinfo/Java7 

DUMP="/home/dwragg/work/Analysis"
PIPE="/home/dwragg/work/Pipeline"
PICARD="/usr/local/bioinfo/src/picard-tools/current"
GATK="/usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.1-1"
REF="/home/dwragg/save/Apis/Apis_mellifera.fa"

# Aim is to take the two high-coverage black bees (18x +) and downsample these
# to around 6x coverage. The downsampled BAMs are then combined to create a
# 12x coverage 'worker' bee. SNP calling is then performed in each of the 6x
# haploid drones, and compared to the results of SNP calling in the 12x diploid
# worker.
SEQAPI="/save/seqapipop/Data/Apis-mellifera/SeqApiPop/BR"

# Soft links to BR20 (18x coerage)
ln -s ${SEQAPI}/BR20-PE_TTAGGC_L008/BR20-PE_TTAGGC_L008_bootstrap.bam \
  ${DUMP}/vcfs/drone-worker/BR20-PE_TTAGGC_L008_bootstrap.bam
ln -s ${SEQAPI}/BR20-PE_TTAGGC_L008/BR20-PE_TTAGGC_L008_bootstrap.bam.bai \
  ${DUMP}/vcfs/drone-worker/BR20-PE_TTAGGC_L008_bootstrap.bam.bai

# Soft links to BR45 (15x coverage)
ln -s ${SEQAPI}/BR45-PE_GATCAG_L008/BR45-PE_GATCAG_L008_bootstrap.bam \
  ${DUMP}/vcfs/drone-worker/BR45-PE_GATCAG_L008_bootstrap.bam
ln -s ${SEQAPI}/BR45-PE_GATCAG_L008/BR45-PE_GATCAG_L008_bootstrap.bam.bai \
  ${DUMP}/vcfs/drone-worker/BR45-PE_GATCAG_L008_bootstrap.bam.bai


# ================================
# Downsample haploids, 5X, 7X, 10X
# ================================
# Needs to be run on both haploid drones
# BR20 probability = 3/18 = 0.17 [result = 6.06x]
# BR45 probability = 3/15 = 0.2 [result = 6.16x]
#SAMPLE="BR20-PE_TTAGGC_L008"
SAMPLE="BR45-PE_GATCAG_L008"
depth="3X"
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/downsample.sh -i ${DUMP}/vcfs/drone-worker/${SAMPLE}_bootstrap.bam \
      -o ${DUMP}/vcfs/drone-worker/${depth}/${SAMPLE}_down.bam \
      -p 0.17
samtools index ${DUMP}/vcfs/drone-worker/${depth}/${SAMPLE}_down.bam


# =================================
# Call SNPs in downsampled haploids
# =================================
# This needs to be run on each of the downsampled haploid samples
#ID="BR20-PE_TTAGGC_L008"
ID="BR45-PE_GATCAG_L008"
depth="7X"
mkdir -p ${DUMP}/vcfs/drone-worker/${depth}/${ID}/vcfs
mkdir -p ${DUMP}/vcfs/drone-worker/${depth}/${ID}/metrics
mkdir -p ${DUMP}/vcfs/drone-worker/${depth}/${ID}/logs
mv ${DUMP}/vcfs/drone-worker/${depth}/${ID}_down.bam ${DUMP}/vcfs/drone-worker/${depth}/${ID}/${ID}_bootstrap.bam
mv ${DUMP}/vcfs/drone-worker/${depth}/${ID}_down.bam.bai ${DUMP}/vcfs/drone-worker/${depth}/${ID}/${ID}_bootstrap.bam.bai
qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
  -o ${DUMP}/logs \
  -e ${DUMP}/logs \
  ${PIPE}/snps.sh -i ${PIPE}/params -p "T" -q "F" -s ${ID} -f ${DUMP}/vcfs/drone-worker/${depth}/${ID}_bootstrap.bam -o ${DUMP}/vcfs/drone-worker/${depth}


# =================================
# Make diploid dones
# =================================
# Merge BAMs (merge both downsampled haploid drones into a diploid drone)
ID="diploid-6X"
depth="3X"
mkdir ${DUMP}/vcfs/drone-worker/${ID}/metrics
mkdir ${DUMP}/vcfs/drone-worker/${ID}/logs
mkdir ${DUMP}/vcfs/drone-worker/${ID}/vcfs
qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 1 \
  -o ${DUMP}/logs \
  -e ${DUMP}/logs \
  ${PIPE}/mergeBAMs.sh -s ${ID} \
    -a ${DUMP}/vcfs/drone-worker/${depth}/BR20-PE_TTAGGC_L008/BR20-PE_TTAGGC_L008_bootstrap.bam  \
    -b ${DUMP}/vcfs/drone-worker/${depth}/BR45-PE_GATCAG_L008/BR45-PE_GATCAG_L008_bootstrap.bam  \
    -o ${DUMP}/vcfs/drone-worker

# Check coverage (run on the diploid sample to confirm coverage correct)
SAMPLE="diploid-6X"
java -d64 -jar ${GATK}/GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R ${REF} \
  -I ${DUMP}/vcfs/drone-worker/${SAMPLE}/${SAMPLE}_bootstrap.bam \
  -o ${DUMP}/vcfs/drone-worker/metrics/${SAMPLE}_GATKcov \
  -ct 2 -ct 5 -ct 8 \
  --omitDepthOutputAtEachBase \
  --omitIntervalStatistics \
  --omitLocusTable \
  -l FATAL


# ================================
# Call SNPs in diploid
# ================================
# This needs to be run on the "diploid done" sample
SAMPLE="diploid-6X"
qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
  -o ${DUMP}/logs \
  -e ${DUMP}/logs \
  ${PIPE}/snps-diploid.sh -i ${PIPE}/params -p "T" -q "F" -s ${SAMPLE} -f ${DUMP}/vcfs/drone-worker/${SAMPLE}/${SAMPLE}_bootstrap.bam -o ${DUMP}/vcfs/drone-worker
# Diploid SNP calling currently in progress for 10X, 14X and 20X
# Coverage needs confirming in both 10X and 20X




# SNPs
# BR20 	= 1059740
# BR45 	= 1079343
# BR 	= 1649963
depth="7X"
diploid="diploid-14X"
gunzip ${DUMP}/vcfs/drone-worker/${depth}/BR20-PE_TTAGGC_L008/vcfs/BR20-PE_TTAGGC_L008_clean.vcf.gz
gunzip ${DUMP}/vcfs/drone-worker/${depth}/BR45-PE_GATCAG_L008/vcfs/BR45-PE_GATCAG_L008_clean.vcf.gz
gunzip ${DUMP}/vcfs/drone-worker/${diploid}/vcfs/${diploid}_clean.vcf.gz



# ==============================================================================
# R
# ==============================================================================

library("VennDiagram")
setwd("/home/dwragg/work/Analysis/vcfs/drone-worker")

depth <- "7X"
diploid <- "diploid-14X"

BR20 <- paste(depth, "/", "/BR20-PE_TTAGGC_L008/vcfs/BR20-PE_TTAGGC_L008_clean.vcf", sep="")
BR45 <- paste(depth, "/", "/BR45-PE_GATCAG_L008/vcfs/BR45-PE_GATCAG_L008_clean.vcf", sep="")
dip <- paste(diploid, "/vcfs/", diploid, "_clean.vcf", sep="")

tmp <- read.table(dip, sep="\t", stringsAsFactors=F)
  varA <- paste(tmp[,1], tmp[,2], sep="-")
  nomA <- rep(diploid, length(varA))
  rm(tmp)
tmp <- read.table(BR20, sep="\t", stringsAsFactors=F)
  varB <- paste(tmp[,1], tmp[,2], sep="-")
  nomB <- rep(paste("BR20-", depth, sep=""), length(varB))
  rm(tmp)
tmp <- read.table(BR45, sep="\t", stringsAsFactors=F)
  varC <- paste(tmp[,1], tmp[,2], sep="-")
  nomC <- rep(paste("BR45-", depth, sep=""), length(varC))
  rm(tmp)
venn.diagram(list(BR=varA, BR20=varB, BR45=varC),
	category.names = c(nomA[1], nomB[1], nomC[1]),
	fill=c("red", "green", "blue"),
	alpha=c(0.5, 0.5, 0.5), cex=0.8,
	cat.fontfamily="sans", cat.fontface=2, cat.cex=0.6,
	fontfamily="sans", fontface=1,
	lty=1, euler.d=TRUE, scaled=TRUE,
	reverse=TRUE, compression='lzw',
	resolution=300, height=5, width=5, units='in', 
	filename=paste(diploid, ".tiff", sep=""));

# ==============================================================================



depth="7X"
diploid="diploid-14X"
bgzip ${DUMP}/vcfs/drone-worker/${depth}/BR20-PE_TTAGGC_L008/vcfs/BR20-PE_TTAGGC_L008_clean.vcf
bgzip ${DUMP}/vcfs/drone-worker/${depth}/BR45-PE_GATCAG_L008/vcfs/BR45-PE_GATCAG_L008_clean.vcf
bgzip ${DUMP}/vcfs/drone-worker/${diploid}/vcfs/${diploid}_clean.vcf













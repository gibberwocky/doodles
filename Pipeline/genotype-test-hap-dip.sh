#!/bin/bash
# IEKNGHER1ujm

# Working directory
DUMP="/home/dwragg/work/Analysis"
PIPE="/home/dwragg/work/Pipeline"

# Ensure vcf files are bgzipped
bgzip -f ${DUMP}/BR20-PE_TTAGGC_L008/vcfs/BR20-PE_TTAGGC_L008_clean.vcf
bgzip -f ${DUMP}/BR45-PE_GATCAG_L008/vcfs/BR45-PE_GATCAG_L008_clean.vcf
bgzip -f ${DUMP}/BR/vcfs/BR_clean.vcf
tabix -fp vcf ${DUMP}/BR20-PE_TTAGGC_L008/vcfs/BR20-PE_TTAGGC_L008_clean.vcf.gz
tabix -fp vcf ${DUMP}/BR45-PE_GATCAG_L008/vcfs/BR45-PE_GATCAG_L008_clean.vcf.gz
tabix -fp vcf ${DUMP}/BR/vcfs/BR_clean.vcf.gz

# Merge sample VCFs into single VCF file
BR=(${DUMP}/B*/vcfs/BR*clean.vcf.gz)
VCFOUT=${DUMP}/tmp/BR-sites.vcf.gz

# Generate VCF file containing merged sites but not genotypes
# Only target chromosome 10
module load bioinfo/bcftools
bcftools merge --regions NC_007079.3 ${BR[@]} \
  | bcftools view --types snps -G -M 2 -O z -o ${VCFOUT} -
tabix -fp vcf ${VCFOUT}

# Generate BAM list
BR_BAM=(${DUMP}/BR/BR*.bam)
ls ${BR_BAM[@]} > ${DUMP}/tmp/BR_BAM.list
# Above file needs checking to ensure it contains correct BAMs

# Genotype samples with GATK using sites from previous step
cd ${DUMP}
qsub -q workq -l mem=4G -l h_vmem=16G -pe parallel_smp 8 -o ${DUMP}/tmp -e ${DUMP}/tmp ${PIPE}/genotype.sh -i ${PIPE}/params -s ${VCFOUT} -b ${DUMP}/tmp/BR_BAM.list -o ${DUMP}/tmp/BR-7079.vcf.gz

# Genotyping using pileup substitutes no coverage with reference allele, not good!
gunzip ${VCFOUT}
sed -i '/^#/ d' ${DUMP}/tmp/BR-sites.vcf
cut -f 1,2 -d '\t' ${DUMP}/tmp/BR-sites.vcf >  ${DUMP}/tmp/BR-sites.bed
module unload bioinfo/bcftools
REF="/home/dwragg/save/Apis/Apis_mellifera.fa"
samtools mpileup -C 50 -q 30 -Q 20 -uDf ${REF} -l ${DUMP}/tmp/BR-sites.bed ${BR_BAM[@]} | bcftools view -cgv - > ${DUMP}/tmp/BR-BCFT.vcf

# Impute
qsub -q workq -l mem=4G -l h_vmem=16G -pe parallel_smp 8 -o ${DUMP}/tmp -e ${DUMP}/tmp ${PIPE}/beagle.sh -i ${DUMP}/tmp/BR-7079.vcf -o ${DUMP}/tmp/BR-7079-beagle.vcf


GATK					BR	BR20	BR45
NC_007079.3:12914239	A	G	0/1	0/1	./.
NC_007079.3:12914294	T	C	1/1	1/1	./.
NC_007079.3:12914320	T	C	1/1	1/1	./.
NC_007079.3:12914325	T	C	1/1	1/1	./.
NC_007079.3:12965464	C	T	0/1	1/1	0/0

PILEUP					BR	BR20	BR45
NC_007079.3:12914239	A	G	0/1	0/1	0/0
NC_007079.3:12914294	T	C	1/1	1/1	1/1
NC_007079.3:12914320	T	C	1/1	1/1	1/1
NC_007079.3:12914325	T	C	1/1	1/1	1/1
NC_007079.3:12965464	C	T	0/1	1/1	0/0

BEAGLE					BR	BR20	BR45
NC_007079.3:12914239	A	G	0|1	1|0	0|1
NC_007079.3:12914294	T	C	1|1	1|1	1|1
NC_007079.3:12914320	T	C	1|1	1|1	1|1
NC_007079.3:12914325	T	C	1|1	1|1	1|1
NC_007079.3:12965464	C	T	0|1	1|1	0|0





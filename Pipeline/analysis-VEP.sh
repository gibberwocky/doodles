#!/bin/bash
# IEKNGHER1ujm

# ! check that VEP and snpEff call same gene annotations!

# =============================================================================
# Ensembl Variant Effect Predictor
# =============================================================================

# SO terms: http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences
# 1) Run for all SNPs to estimate Pn/Ps per gene per population
# 2) Run for non-reference private SNPs per population
# 3) Run for signature of selection results


# Load modules
module load bioinfo/bcftools

# Paths to main folders
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PERL=${PIPE}/Perl
VCFT="/usr/local/bioinfo/src/vcftools/vcftools_0.1.12a/bin"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"

# Pre-imputation haploid JFM OTH data
SITES=${DUMP}/vcfs/HARP-JFM-OTH-15012015
VCF="JFMOTH_preImp_chrs16_qc_beagle"

# =============================================================================
# Generate INFO VCF files for each subset of samples per population
# The AC and AN counts refer to alt allele and total alleles within the
# population. The AF field however relates to the whole data before sample
# sub-setting. The frequencies will be used in combination with the VEP
# output to calculate in R statistics per gene relating to pN/pS, maybe.
# =============================================================================
POPS=('JFM' 'OTH')
for POP in ${POPS[@]}
do
 # Output population-specific INFO (AC and AN only; AF relates to whole data)
  bcftools view -H -G -O v -o ${SITES}/${VCF}_${POP}_info.vcf \
    -S ${DUMP}/vcfs/lists/${POP}_samples_fixed.list \
    ${SITES}/${VCF}.vcf.gz
  # Reduce down to just chromosome annotation and re-calculated AF (AC/AN)
  sed 's/;/\t/g' ${SITES}/${VCF}_${POP}_info.vcf | cut -f1-5,11-12 - | sed 's/AC=//g' - | sed 's/AN=//g' - | awk -F $'\t' 'BEGIN {OFS = FS} {print $1,  $2, $3, $4, $5, $6/$7}' - > ${SITES}/${VCF}_${POP}_AF.vcf
  # Remove surplous file
  ${SITES}/${VCF}_${POP}_info.vcf
done


# Relate to above but not to be run at same time, this parses a GTF file
# and adds column with gene id
GTF=/save/dwragg/Apis/Apis_mellifera.GCA_000002195.1.25.gtf
awk -F'\t' 'BEGIN {OFS = FS}; {split($9, g, ";"); printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $7, g[1]);}' ${GTF} | sed 's/gene_id//g' - | sed 's/\"//g' - | sed 's/\t */\t/g' - > ${GTF}.clean




mkdir -p ${DUMP}/logs/VEP
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
  -o ${DUMP}/logs/VEP \
  -e ${DUMP}/logs/VEP \
  ${PIPE}/vep.sh \
    -i ${SITES}/${VCF}_${POP}_info.vcf  \
    -o ${DUMP}/VEP/${VCF}.vcf 


sed '/^##/ d' ${DUMP}/VEP/${VCF}.vcf_vep_output.txt | cut -f1-2,4,7 - > ${DUMP}/VEP/${VCF}_consequences.vep

# =============================================================================
# Private non-reference SNPs
# =============================================================================
POPS=('JFM' 'OTH')
for POP in ${POPS[@]}
do
  bcftools view -M 2 -O v --private \
    -S ${DUMP}/vcfs/lists/${POP}_samples_fixed.list \
    -o ${SITES}/${VCF}_${POP}_private.vcf \
    ${SITES}/${VCF}.vcf.gz
  mkdir -p ${DUMP}/logs/VEP-${POP}-private
  qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/VEP-${POP}-private \
    -e ${DUMP}/logs/VEP-${POP}-private \
    ${PIPE}/vep.sh \
      -i ${SITES}/${VCF}_${POP}_private.vcf  \
      -o ${DUMP}/VEP/${VCF}_${POP}_private.vcf 
done



# =============================================================================
# All non-reference SNPs
# =============================================================================
POPS=('JFM' 'OTH')
for POP in ${POPS[@]}
do
  bcftools view -M 2 -O v --min-ac 3 \
    -S ${DUMP}/vcfs/lists/${POP}_samples_fixed.list \
    -o ${SITES}/${VCF}_${POP}_nonRef.vcf \
    ${SITES}/${VCF}.vcf.gz
  mkdir -p ${DUMP}/logs/VEP-${POP}-nonRef
  qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/VEP-${POP}-nonRef \
    -e ${DUMP}/logs/VEP-${POP}-nonRef \
    ${PIPE}/vep.sh \
      -i ${SITES}/${VCF}_${POP}_nonRef.vcf  \
      -o ${DUMP}/VEP/${VCF}_${POP}_nonRef.vcf 
done










# =============================================================================
# 3. SNPs found to be statistically significant (ZHp)
# =============================================================================
POP="JFM"
hp="pos"
HITS=${DUMP}/selection/hp/${POP}-${hp}-merged-multi.bed
SITES=${DUMP}/vcfs/HARP-JFM-OTH-07112014

# Change BED chromosome names to NCBI accessions
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/bedchr-to-accession.sh \
      -i ${HITS} 

# Retrieve SNPs
bcftools view -M 2 -O v \
  -S ${DUMP}/vcfs/lists/${POP}_samples.list \
  -o ${DUMP}/VEP/${POP}-hp-${hp}.vcf \
  -R ${HITS} \
  ${SITES}/S2_hapJFMOTH_beagle.vcf.gz

# Flip chromosomes back to numbers for ensembl
sed -i '/^#CHROM/s/_/-/g' ${DUMP}/VEP/${POP}-hp-${hp}.vcf
${PIPE}/accession-to-chr.sh -i ${DUMP}/VEP/${POP}-hp-${hp}.vcf

# Run VEP
mkdir -p ${DUMP}/logs/VEP-${POP}
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/VEP-${POP} \
    -e ${DUMP}/logs/VEP-${POP} \
    ${PIPE}/vep.sh \
      -i ${DUMP}/VEP/${POP}-hp-${hp}.vcf \
      -o ${DUMP}/VEP/${POP}-hp-${hp}

# =============================================================================
# 3. SNPs found to be statistically significant (Fst)
# =============================================================================
HITS=${DUMP}/selection/vcft-fst/fst.bed
SITES=${DUMP}/vcfs/HARP-JFM-OTH-07112014

# Merge overlapping features
bedtools merge -c 1 -o count -i ${HITS} > ${HITS}.merged

# Change BED chromosome names to NCBI accessions
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/bedchr-to-accession.sh \
      -i ${HITS}.merged

# Retrieve SNPs
bcftools view -M 2 -O v \
  -o ${DUMP}/VEP/JFMOTH-fst.vcf \
  -R ${HITS}.merged \
  ${SITES}/S2_hapJFMOTH_beagle.vcf.gz

# Flip chromosomes back to numbers for ensembl
sed -i '/^#CHROM/s/_/-/g' ${DUMP}/VEP/JFMOTH-fst.vcf
${PIPE}/accession-to-chr.sh -i ${DUMP}/VEP/JFMOTH-fst.vcf

# Run VEP
mkdir -p ${DUMP}/logs/VEP-${POP}
qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/VEP-${POP} \
    -e ${DUMP}/logs/VEP-${POP} \
    ${PIPE}/vep.sh \
      -i ${DUMP}/VEP/JFMOTH-fst.vcf \
      -o ${DUMP}/VEP/JFMOTH-fst






# =============================================================================
# 4. Feature intersects between JFM and OTH Hp results
# =============================================================================


cat ${DUMP}/selection/hp/JFM-neg-merged-multi.bed \
  ${DUMP}/selection/hp/JFM-pos-merged-multi.bed \
  > ${DUMP}/selection/hp/JFM-hp.bed
cat ${DUMP}/selection/hp/OTH-neg-merged-multi.bed \
  ${DUMP}/selection/hp/OTH-pos-merged-multi.bed \
  > ${DUMP}/selection/hp/OTH-hp.bed

bedtools intersect -loj \
  -b ${DUMP}/selection/hp/JFM-hp.bed -a ${DUMP}/selection/hp/OTH-hp.bed 
























#!/bin/bash
# IEKNGHER1ujm

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


# ==============================================================================
# Fst
# ==============================================================================

# Copy sample lists
cp ${DUMP}/vcfs/lists/JFM_samples.list ${DUMP}/selection/vcft-fst/JFM_samples.list
sed -i 's/_/-/g' ${DUMP}/selection/vcft-fst/JFM_samples.list
cp ${DUMP}/vcfs/lists/OTH_samples.list ${DUMP}/selection/vcft-fst/OTH_samples.list
sed -i 's/_/-/g' ${DUMP}/selection/vcft-fst/OTH_samples.list

# Calculate windowed Fst (Weir & Cockerham 1984)
${VCFT}/vcftools --gzvcf ${SITES}/${VCF}.vcf.gz \
  --weir-fst-pop ${DUMP}/selection/vcft-fst/JFM_samples.list \
  --weir-fst-pop ${DUMP}/selection/vcft-fst/OTH_samples.list \
  --fst-window-size 50000 \
  --fst-window-step 10000 \
  --out ${DUMP}/selection/vcft-fst/${VCF}
# 3159533-26


# ==============================================================================
# Hp 
# ==============================================================================

# Copy sample lists
cp ${DUMP}/vcfs/lists/JFM_samples.list ${DUMP}/selection/hp/JFM_samples.list
sed -i 's/_/-/g' ${DUMP}/selection/hp/JFM_samples.list
cp ${DUMP}/vcfs/lists/OTH_samples.list ${DUMP}/selection/hp/OTH_samples.list
sed -i 's/_/-/g' ${DUMP}/selection/hp/OTH_samples.list

POP="OTH"
bcftools view -S ${DUMP}/selection/hp/${POP}_samples.list \
  -O v -o ${DUMP}/selection/hp/${VCF}-${POP}.vcf \
  ${SITES}/${VCF}.vcf.gz

# Fill AN if missing:
perl ${VCFT}/fill-an-ac ${DUMP}/selection/hp/${VCF}-${POP}.vcf \
  > ${DUMP}/selection/hp/${VCF}-${POP}-AN.vcf
mv -f ${DUMP}/selection/hp/${VCF}-${POP}-AN.vcf \
  ${DUMP}/selection/hp/${VCF}-${POP}.vcf

# Windowed Hp analysis (VCF CHR WIN_SIZE WIN_STEP)
perl ${PERL}/VCF-Hp.pl ${DUMP}/selection/hp/${VCF}-${POP}.vcf 16 50000 10000
perl ${PERL}/VCF-Hp.pl ${DUMP}/selection/hp/${VCF}.vcf 16 50000 10000

# Output average Hp
awk 'NR == 1 { sum=0 } { sum+=$10;} END {printf "Average: %f\n", sum/NR}' \
    ${DUMP}/selection/hp/${VCF}-${POP}.vcf.hp




# For a combined JFM+OTH analysis just gunzip the main VCF
bcftools view -O v -o ${DUMP}/selection/hp/${VCF}.vcf \
  ${SITES}/${VCF}.vcf.gz
perl ${VCFT}/fill-an-ac ${DUMP}/selection/hp/${VCF}.vcf \
  > ${DUMP}/selection/hp/${VCF}-AN.vcf
mv -f ${DUMP}/selection/hp/${VCF}-AN.vcf \
  ${DUMP}/selection/hp/${VCF}.vcf
perl ${PERL}/VCF-Hp.pl ${DUMP}/selection/hp/${VCF}.vcf 16 5000 1000






# Check for any intersecting regions after processing results in R
OUT=${DUMP}/selection/results/raw
cut -f1-3 ${OUT}/OTH-pos.bed.clean > ${OUT}/OTH-pos.tmp
cut -f1-3 ${OUT}/OTH-neg.bed.clean > ${OUT}/OTH-neg.tmp
cat ${OUT}/OTH-neg.tmp ${OUT}/OTH-neg.tmp > ${OUT}/OTH.tmp.bed
cut -f1-3 ${OUT}/JFM-pos.bed.clean > ${OUT}/JFM-pos.tmp
cut -f1-3 ${OUT}/JFM-neg.bed.clean > ${OUT}/JFM-neg.tmp
cat ${OUT}/JFM-neg.tmp ${OUT}/JFM-neg.tmp > ${OUT}/JFM.tmp.bed
bedtools intersect -a ${OUT}/JFM.tmp.bed -b ${OUT}/OTH.tmp.bed


# ==============================================================================
# ADMIXTURE / fastSTRUCTURE
# ==============================================================================

SITES=${DUMP}/vcfs/HARP-JFM-OTH-15012015
VCF="Harpur-dipJFMOTH_preImp_chrs16_qc"
#VCF="Harpur-hapJFMOTH_preImp_chrs16_qc"
#VCF="Harpur-hapJFMOTH_preImp_chrs16_qc_beagle"
mkdir -p ${DUMP}/vcfs/plink/${VCF}

# Make BED file (subset of samples)
#${PLINK}/plink --vcf ${SITES}/${VCF}.vcf \
#  --keep-fam ${SITES}/Harpur-JFMOTH-Liu.fam \
#  --keep-allele-order \
#  --a2-allele ${SITES}/${VCF}.vcf 4 3 '#' \
#  --allow-no-sex \
#  --allow-extra-chr \
#  --mind 1 \
#  --geno 0.1 \
#  --maf 0.05 \
#  --chr-set 16 \
#  --chr 1-16 \
#  --set-missing-snp-ids @:#[Amel4-5] \
#  --out ${DUMP}/vcfs/plink/${VCF}/${VCF}-no-Wallberg \
#  --make-bed

# Make BED file
${PLINK}/plink --vcf ${SITES}/${VCF}.vcf \
  --keep-allele-order \
  --a2-allele ${SITES}/${VCF}.vcf 4 3 '#' \
  --allow-no-sex \
  --allow-extra-chr \
  --mind 1 \
  --geno 0.1 \
  --maf 0.05 \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/plink/${VCF}/${VCF} \
  --make-bed

# Before running ADMIXTURE first remove SNPs in high LD
#${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${VCF}/${VCF} \
#  --keep-allele-order \
#  --allow-no-sex \
#  --allow-extra-chr \
#  --mind 1 \
#  --geno 0.1 \
#  --maf 0.01 \
#  --chr-set 16 \
#  --chr 1-16 \
#  --set-missing-snp-ids @:#[Amel4-5] \
#  --out ${DUMP}/vcfs/plink/${VCF}/${VCF} \
#  --indep-pairwise 50 10 0.1
#${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${VCF}/${VCF} \
#  --keep-allele-order \
#  --allow-no-sex \
#  --allow-extra-chr \
#  --mind 1 \
#  --geno 0.1 \
#  --maf 0.01 \
#  --bp-space 5000 \
#  --chr-set 16 \
#  --chr 1-16 \
#  --set-missing-snp-ids @:#[Amel4-5] \
#  --extract ${DUMP}/vcfs/plink/${VCF}/${VCF}.prune.in \
#  --out ${DUMP}/vcfs/plink/${VCF}/${VCF}-prune \
#  --make-bed

# Harpur-dipJFMOTH_preImp_chrs16_qc-prune
# 32618 variants and 69 samples pass filters and QC

# Run Admixture on cluster
k=6 # Min K
K=6 # Max K
N=10 # Number of bootstraps
SUFFIX=""
mkdir -p ${DUMP}/vcfs/admixture/${VCF}${SUFFIX}
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/admixture.sh -i ${DUMP}/vcfs/plink/${VCF}/${VCF}${SUFFIX} \
      -k ${k} -K ${K} -n ${N} -o ${DUMP}/vcfs/admixture/${VCF}${SUFFIX}/ \
      -s ${VCF}${SUFFIX}
# Copy fam file to same folder incase needed for plotting
cp ${DUMP}/vcfs/plink/${VCF}/${VCF}${SUFFIX}.fam \
  ${DUMP}/vcfs/admixture/${VCF}${SUFFIX}
# Identify optimal K
grep -h CV ${DUMP}/vcfs/admixture/${VCF}${SUFFIX}/${VCF}-*.out

# Run fastSTRUCTURE on cluster
k=1 # Min K
K=6 # Max K
CV=0
mkdir -p ${DUMP}/vcfs/fastSTRUCTURE/${VCF}${SUFFIX}
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs -e ${DUMP}/logs  \
    ${PIPE}/fastStructure.sh -i ${DUMP}/vcfs/plink/${VCF}/${VCF}${SUFFIX} \
      -k ${k} -K ${K} -n ${CV} \
      -o ${DUMP}/vcfs/fastSTRUCTURE/${VCF}${SUFFIX}/${VCF}${SUFFIX}

#chooseK.py --input=${DUMP}/vcfs/fastStructure/${VCF}${SUFFIX}/${VCF}${SUFFIX}
#distruct.py -K 6 --input=${OUT}/fastStructure/test \
#  --popfile=${OUT}/popfile \
#  --output=${OUT}/fastStructure/test-K6.svg


# =============================================================================
# Plink to R - for phylogenetic analysis
# =============================================================================

# This assumes that the Plink file used in the ADMIXTURE analysis is intacted
# and in the ${DUMP}/vcfs/plink folder, and that the ADMIXTURE results have
# been transferred to their own folder in ${DUMP}/vcfs/admixture, including
# the plink fam file.

SITES=${DUMP}/vcfs/HARP-JFM-OTH-15012015
VCF="Harpur-JFMOTH_chrs16_qc_beagle"
INDIR=${DUMP}/vcfs/plink/${VCF}

# FAM to PHENO (id	sex	 phene)
cut -d' ' -f1,3,4 ${INDIR}/${VCF}.fam > ${INDIR}/${VCF}-R.tmp
header=('id' 'sex' 'pop')
echo ${header[@]} | cat - ${INDIR}/${VCF}-R.tmp > ${INDIR}/${VCF}-R.pheno
rm ${INDIR}/${VCF}-R.tmp

# Format Plink for GenABEL import
${PLINK}/plink --bfile ${DUMP}/vcfs/plink/${VCF}/${VCF} \
  --keep-allele-order \
  --allow-no-sex \
  --allow-extra-chr \
  --mind 1 \
  --geno 0.1 \
  --maf 0.05 \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --recode 12 \
  --out ${DUMP}/vcfs/plink/${VCF}/${VCF}-R

# Remove genetic distance from map file
cut -f1,2,4 ${DUMP}/vcfs/plink/${VCF}/${VCF}-R.map \
  > ${DUMP}/vcfs/plink/${VCF}/${VCF}-R.tmp
mv ${DUMP}/vcfs/plink/${VCF}/${VCF}-R.tmp ${DUMP}/vcfs/plink/${VCF}/${VCF}-R.map

# Proceed in R









# =============================================================================
# Haplotype analysis
# =============================================================================

# require separate file per chromosome
SITES=${DUMP}/vcfs/HARP-JFM-OTH-15012015
VCF="JFMOTH_preImp_chrs16_qc_beagle"

# Homozygous SNPs only
#bcftools view -g ^het -O z -o ${SITES}/${VCF}_noHet.vcf.gz \
#  ${SITES}/${VCF}.vcf.gz
#tabix -fp ${SITES}/${VCF}_noHet.vcf.gz


mkdkr -p ${DUMP}/vcfs/haploblocks/CHR
POP="OTH"
CHRS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16')
mkdir -p ${DUMP}/vcfs/haploblocks/${POP}
for CHR in ${CHRS[@]}
do
  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/vcfs/haploblocks/logs \
    -e ${DUMP}/vcfs/haploblocks/logs \
    ${PIPE}/haps.sh -v ${SITES}/${VCF}.vcf.gz \
      -o ${DUMP}/vcfs/haploblocks \
      -p ${POP} -c ${CHR} -s ${DUMP}/selection/hp/${POP}_samples.list
done


# Change allele states:  0 = missing, 1 = reference, 2 = alternative
#FILES=(${DUMP}/vcfs/haploblocks/*/*.hap)
#for ID in ${FILES[@]}
#do
#  sed 's/1/2/g' ${ID} > ${ID}.012
#  sed 's/0/1/g' ${ID} > ${ID}.012
#  sed 's/\./0/g' ${ID} > ${ID}.012
#done



# =============================================================================
# iHs / XP-EHH - data must be phased and have no missing genotypes
# Use Beagle to impute missing genotypes on JFMOTH_preImp_chrs16_qc
# Genetic position based on Apis mellifera recombination = ~ 37 cM/Mb
# =============================================================================

# -----------------------------------------------------------------------------
# XP-EHH
# -----------------------------------------------------------------------------

POP1="JFM"
POP2="OTH"
CHRS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16')
for CHR in ${CHRS[@]}
do
  qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/selection/selscan/xpehh/logs \
    -e ${DUMP}/selection/selscan/xpehh/logs \
    ${PIPE}/selscan.sh -w 100000 -g 200000 -s 20000 -x 1000000 \
      -h ${DUMP}/vcfs/haploblocks/${POP1}/${POP1}-${CHR}.hap \
      -r ${DUMP}/vcfs/haploblocks/${POP2}/${POP2}-${CHR}.hap \
      -m ${DUMP}/vcfs/haploblocks/CHR/${CHR}.map \
      -o ${DUMP}/selection/selscan/xpehh/${POP1}-${POP2}-${CHR}
done

# XP-EHH post-processing (do not repeat more than once!!!!)
for CHR in ${CHRS[@]}
do
  # Remove header line
  sed -i 1d ${DUMP}/selection/selscan/xpehh/${POP1}-${POP2}-${CHR}.xpehh.out
  # Add chromosome field
  sed -i "s/^/${CHR}\t/" \
    ${DUMP}/selection/selscan/xpehh/${POP1}-${POP2}-${CHR}.xpehh.out
done
cat ${POP1}-${POP2}-*.xpehh.out > ${POP1}-${POP2}.xpehh

# -----------------------------------------------------------------------------
# iHS
# -----------------------------------------------------------------------------

POP="JFM"
for CHR in ${CHRS[@]}
do
  qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/selection/selscan/ihs/logs \
    -e ${DUMP}/selection/selscan/ihs/logs \
    ${PIPE}/ihs.sh \
      -h ${DUMP}/vcfs/haploblocks/${POP}/${POP}-${CHR}.hap \
      -m ${DUMP}/vcfs/haploblocks/CHR/${CHR}.map \
      -o ${DUMP}/selection/selscan/ihs/${POP}-${CHR}
done

# iHS post-processing
POP="OTH"
for CHR in ${CHRS[@]}
do
  # Normalize in windows as per docs
  norm --files ${DUMP}/selection/selscan/ihs/${POP}-${CHR}.ihs.out
  # Add chromosome field
  sed -i "s/^/${CHR}\t/" \
    ${DUMP}/selection/selscan/ihs/${POP}-${CHR}.ihs.out.100bins.norm
done
cat ${DUMP}/selection/selscan/ihs/${POP}-*.ihs.out.100bins.norm \
    > ${DUMP}/selection/selscan/ihs/${POP}.ihs






# =============================================================================
# Processing of 'significant.bed' results from selection analyses
# This is integrated into the Plotting script so may be redundant
# =============================================================================

BED2="/home/dave/Downloads/bedtools2/bin"
OUT="/home/dave/Copy/HoneyBee/Analyses/selection/Results"

# Sort and merge the bed file
${BED2}/bedtools sort -i ${OUT}/significant.bed \
  > ${OUT}/significant.srt.bed
${BED2}/bedtools merge -i ${OUT}/significant.srt.bed \
  -d 5000 -c 1,4,5 -o count,absmax,distinct \
  | sort -k1,1n -k2,2n - > ${OUT}/significant.merged.bed






# =============================================================================
# Check nucleotide diversity
# =============================================================================

SITES=${DUMP}/vcfs/HARP-JFM-OTH-15012015
VCF="Harpur-hapJFMOTH_preImp_chrs16_qc.vcf.gz"

POP="OTH"

# Calculate windowed Fst (Weir & Cockerham 1984)
${VCFT}/vcftools --gzvcf ${SITES}/${VCF}.vcf.gz \
  --keep ${DUMP}/selection/hp/${POP}_samples.list \
  --window-pi 50000 \
  --window-pi-step 10000 \
  --out ${DUMP}/selection/vcft-nuc/${VCF}-${POP}

# Output average
awk 'NR == 1 { sum=0 } { sum+=$5;} END {printf "Average: %f\n", sum/NR}' \
    ${DUMP}/selection/vcft-nuc/${VCF}-${POP}*.pi
#JFM = 0.003845
#OTH = 0.003482




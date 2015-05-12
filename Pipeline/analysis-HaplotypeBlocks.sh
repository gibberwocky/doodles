#!/bin/bash
# IEKNGHER1ujm

# Set-up
module load bioinfo/bcftools
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
VCFT="/usr/local/bioinfo/src/vcftools/vcftools_0.1.12a/bin"


# Use input file from Hp analysis
#VCF="S2_hapJFMOTH_beagle.vcf.gz"
VCF="S2_hapJFMOTH_preImp.vcf.gz"
cp ${DUMP}/vcfs/HARP-JFM-OTH-07112014/${VCF} ${DUMP}/vcfs/haploblocks/${VCF}
tabix -fp vcf ${DUMP}/vcfs/haploblocks/${VCF}

# Extract chromosomes 1 to 16
#FILE="JFM-OTH"
FILE="JFM-OTH-preImp"
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3' \
  -o ${DUMP}/vcfs/haploblocks/${FILE}.vcf \
  ${DUMP}/vcfs/haploblocks/${VCF}

# Change from accession to chromosome naming
${PIPE}/accession-to-chr.sh -i ${DUMP}/vcfs/haploblocks/${FILE}.vcf
sed -i '/^#CHROM/s/_/-/g' ${DUMP}/vcfs/haploblocks/${FILE}.vcf

# Make BED file
${PLINK}/plink --vcf ${DUMP}/vcfs/haploblocks/${FILE}.vcf \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr 1-16 \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/haploblocks/${FILE} \
  --make-bed

# Create family files to filter on in Plink
head --l 30 ${DUMP}/vcfs/haploblocks/${FILE}.fam > ${DUMP}/vcfs/haploblocks/JFM.fam
tail --l 30 ${DUMP}/vcfs/haploblocks/${FILE}.fam > ${DUMP}/vcfs/haploblocks/OTH.fam

# Calculated haplotype blocks for nominated chromosome (det output file)
# Perform once without sample filtering, and once for each population
POP="OTH"
CHROMOSOME=6
qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
  -o ${DUMP}/vcfs/haploblocks \
  -e ${DUMP}/vcfs/haploblocks \
  ${PIPE}/haploblocks.sh -c ${CHROMOSOME} \
    -b ${DUMP}/vcfs/haploblocks/${FILE} \
    -o ${DUMP}/vcfs/haploblocks/${POP}-${CHROMOSOME} \
    -f ${DUMP}/vcfs/haploblocks/${POP}.fam


# Copy the *.det files to laptop for plotting (well, that's the idea)


# To iIsolate blocks in an area of interest:
chr=6
pos1=6205001
pos2=6226000
awk -v VAR1=${chr} -v VAR2=${pos1} -v VAR3=${pos2} \
  '$1 ==VAR1 && $2 >=VAR2 && $2 <=VAR3' \
  ${DUMP}/vcfs/haploblocks/${FILE}.blocks.det \
  > ${DUMP}/vcfs/haploblocks/${FILE}.blocks.tmp
head --l 1 ${DUMP}/vcfs/haploblocks/${FILE}.blocks.det | cat - ${DUMP}/vcfs/haploblocks/${FILE}.blocks.tmp > ${DUMP}/vcfs/haploblocks/${FILE}.blocks.hit
rm *.tmp


# Calculated runs of homozygosity
CHROMOSOME=6
${PLINK}/plink --bfile ${DUMP}/vcfs/haploblocks/${FILE} \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --chr ${CHROMOSOME} \
  --set-missing-snp-ids @:#[Amel4-5] \
  --out ${DUMP}/vcfs/haploblocks/${FILE} \
  --homozyg group extend
# Copy fam file to same folder incase needed for plotti






# Calculate haplotype blocks (hap output file)
# s is the start column in the VCF file of the population
# n is the number of samples (cols) for the population
# x (T/F) indicates whether to include sites with missing genotypes, defualt=T
# h (T/F) indicates whether to include sites with heterozygous genotypes, defualt=F
# script will not overwrite CHR files if the exist already

# *** preImp needs renaming as its subhet ***

VCF="S2_hapJFMOTH_beagle.vcf.gz"
VCF="S2_hapJFMOTH_preImp.vcf.gz"
VCF="S2_hapJFMOTH_beagle_subhet.vcf.gz"
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/vcfs/haploblocks/haps/subhet \
    -e ${DUMP}/vcfs/haploblocks/haps/subhet \
    ${PIPE}/haps.sh -v ${DUMP}/vcfs/haploblocks/${VCF} \
      -o ${DUMP}/vcfs/haploblocks/haps/subhet \
      -p OTH -s 40 -n 30 -h F

# Change allele states:  0 = missing, 1 = reference, 2 = alternative
FILES=(${DUMP}/vcfs/haploblocks/haps/subhet/*/*.hap)
for ID in ${FILES[@]}
do
  sed -i 's/1/2/g' ${ID}
  sed -i 's/0/1/g' ${ID}
  sed -i 's/\./0/g' ${ID}
done


# Replace heterozygous values with missing in beagle VCF file
sed -i 's/0|1/.|./g' ${DUMP}/vcfs/haploblocks/S2_hapJFMOTH_beagle_subhet.vcf
sed -i 's/1|0/.|./g' ${DUMP}/vcfs/haploblocks/S2_hapJFMOTH_beagle_subhet.vcf  
bgzip -f ${DUMP}/vcfs/haploblocks/S2_hapJFMOTH_beagle_subhet.vcf
tabix -fp vcf ${DUMP}/vcfs/haploblocks/S2_hapJFMOTH_beagle_subhet.vcf.gz 




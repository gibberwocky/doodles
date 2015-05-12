#!bin/bash

# Modify environment for Java
module load bioinfo/Java7 
module load bioinfo/bcftools


SNPEFF="/usr/local/bioinfo/src/SnpEff/snpEff_3.6c"
DUMP="/home/dwragg/work/Analysis"
PIPE="/home/dwragg/work/Pipeline"

# To run SnpEff on all SNPs identified
# Strip out unplaced contigs
bcftools view --types snps -M 2 -O v \
  -o ${DUMP}/vcfs/EurBee-NCBI.vcf \
  -r NC_007070.3,\
NC_007071.3,\
NC_007072.3,\
NC_007073.3,\
NC_007074.3,\
NC_007075.3,\
NC_007076.3,\
NC_007077.3,\
NC_007078.3,\
NC_007079.3,\
NC_007080.3,\
NC_007081.3,\
NC_007082.3,\
NC_007083.3,\
NC_007084.3,\
NC_007085.3,\
NC_001566.1,\
  ${DUMP}/vcfs/EurBee.vcf.gz
# Then run SnpEff
java -jar ${SNPEFF}/snpEff.jar -v amel4.5 ${DUMP}/vcfs/EurBee-NCBI.vcf > ${DUMP}/vcfs/EurBee-NCBI-eff.vcf





# ==============================================================================
# TO ANALYZE A SUBSET OF DATA, EG FROM A SELECTIVE SWEEP ANALYSIS
# ==============================================================================

# First create a VCF file filtered on the results from the Hp test
# To do this requires writing a list of candidate regions to file, to then
# use as a filter on the VCF file 
# Can export table of hits from R, then merge overlapping features with
# bedtools:


bedtools merge -d 10000 -i ${DUMP}/vcfs/Fst-10k-25step.out.bed > ${DUMP}/vcfs/Fst.bed


# Then extract SNPs in regions:


bcftools view --types snps -M 2 -O v \
  -o ${DUMP}/vcfs/EurBee-Fst.vcf \
  -R ${DUMP}/vcfs/Fst.bed \
  ${DUMP}/vcfs/EurBee.vcf.gz


# Then run SnpEff:

java -jar ${SNPEFF}/snpEff.jar -v amel4.5 ${DUMP}/vcfs/EurBee-Fst.vcf > ${DUMP}/vcfs/EurBee-Fst-eff.vcf


# Just genes (remove chromosome IDs, duplicates, and sort)
java -jar ${SNPEFF}/SnpSift.jar extractFields \
  ${DUMP}/vcfs/EurBee-Fst-eff.vcf \
  EFF[*].GENE \
  | sed 's/\t/\n/g' - \
  | sed '/^#/d' - \
  | sed 's/^[^:]*://' - | sort - | uniq - \
  > Fst-Hits.genes


# Delete all lines except those matching LG11 ID, requires previous tool to not
# strip out the chromosome IDs
sed '/^[^NC_007080.3]*NC_007080.3/!d' Hits-1-Hp.genes.tmp > Hits-1-Hp.LG11.genes

# Compare two sets of genes, output genes unique to file 1
comm Hits-1-Hp.genes Hits-2-Hp.genes \
  | awk -F "\t" '{print $1}' - \
  | sed '/^$/d' - \
  > Hits-1-Hp.uniq.genes

# Compare two sets of genes, output genes unique to file 2
comm Hits-1-Hp.genes Hits-2-Hp.genes \
  | awk -F "\t" '{print $2}' - \
  | sed '/^$/d' - \
  > Hits-2-Hp.uniq.genes





# Subset VCF by particular variants (not useful at present)
bcftools view --types snps -M 2 -O v \
  -o ${DUMP}/vcfs/Hits-1-Hp-nonSyn.vcf \
  -i 'INFO/EFF[*]~"NON_SYN"' \
  Hits-1-Hp-eff.vcf

# Tabulate effect details per SNP (not useful at present)
java -jar ${SNPEFF}/SnpSift.jar extractFields -s "," \
  ${DUMP}/vcfs/Hits-1-Hp-eff.vcf \
  CHROM POS EFF[*].EFFECT EFF[*].IMPACT EFF[*].GENE GEN[*].GT \
  | head




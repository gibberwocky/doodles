#!bin/bash

# 1) Extract gene region from VCF to new VCF file
# 2) Upload VCF to Ensembl VEP
# 3) Download results (CSV) and copy into Calc
# 4) Filter on gene and missense SNPs
# 5) Enter protein sequence and SNP details into SIFT

module load bioinfo/bcftools
bcftools view -r NC_007070.3:18531763-18535982 --types snps -G -M 2 -O v -o test-1.vcf EurBee.vcf.gz
sed -i 's/NC_007070.3/1/g' test-1.vcf


>Lgr3
NC_007070.3:18514554-18520902
>Or160
NC_007070.3:18522235-18525065
>HRG
NC_007070.3:18531763-18535982






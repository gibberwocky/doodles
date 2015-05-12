#!/bin/bash
# IEKNGHER1ujm

# Core Information


# Paths to files on seqapipop and for creating links:
paths.sh

# Details on generating a set of master sites:
sites.sh

# Genotyping and imputing in a population:
sites.sh
	(depends imp-geno.sh)
	(diploid-drone.sh [see notes])



# Admixture/ Population analysis
Requires running on diploid dataset, SNPs common to AMCO and SeqApiPop
Use same dataset with Admixture and phylogenetics

# Signatures of selection
Only SNPs detected within haploid dataset relevant here
Use seqapipop master sites (JFMOTH_raw.vcf.gz) to subset S2_hapJFMOTH_beagle.vcf 
This to ensure imputed genotypes are included 


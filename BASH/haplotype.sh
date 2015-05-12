#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# Set-up
module load bioinfo/bcftools
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"

#VCF=${DUMP}/selection/ehh/S2_hapJFMOTH_beagle.vcf.gz
VCF=${DUMP}/selection/ehh/S2_hapJFMOTH_preImp.vcf.gz
FILE=${DUMP}/selection/ehh/JFM-OTH.vcf

CHRS=('NC_007070.3' 'NC_007071.3' 'NC_007072.3' 'NC_007073.3' \
  'NC_007074.3' 'NC_007075.3' 'NC_007076.3' 'NC_007077.3' 'NC_007078.3' \
  'NC_007079.3' 'NC_007080.3' 'NC_007081.3' 'NC_007082.3' 'NC_007083.3' \
  'NC_007084.3' 'NC_007085.3' 'NC_001566.1')

# Chromosomal loop
function join { local IFS="$1"; shift; echo "$*"; }

for CHR in ${CHRS[@]}
do

	# Extract chromosomes, exclude sites with uncalled genotype, drop header
	bcftools view -H -g ^het -g ^miss -O v \
	  --regions ${CHR} -o ${FILE} ${VCF}

	# Extract header and save sample IDs
	bcftools view -h  -O v \
	  -o ${DUMP}/selection/ehh/header.vcf \
	  ${VCF}
	tail --l 1 ${DUMP}/selection/ehh/header.vcf \
	  > ${DUMP}/selection/ehh/IDs
	rm ${DUMP}/selection/ehh/header.vcf

	# Number of samples
	POP="OTH"
	START=40
	N=30
	mkdir -p ${DUMP}/selection/ehh/${POP}
	HAPFILE=${DUMP}/selection/ehh/${POP}/${POP}-${CHR}.hap

	rm -f ${HAPFILE}
	touch ${HAPFILE}

	for((I=$((${START}));I<=$((${START}+${N}-1));I++));
	do 
	  # Cut n paste
	  cut -f${I} ${FILE} > ${DUMP}/selection/ehh/test.${I}
	  cut -c 1 ${DUMP}/selection/ehh/test.${I} > ${DUMP}/selection/ehh/test.${I}a
	  readarray tmp < ${DUMP}/selection/ehh/test.${I}a
	  # Append data
	  echo ${tmp[@]} >> ${HAPFILE}
	  # Clean up
	  rm ${DUMP}/selection/ehh/test.${I}
	  rm ${DUMP}/selection/ehh/test.${I}a
	done


	# Extract mapfile details
	mkdir -p ${DUMP}/selection/ehh/CHR
	cut -f 1-3 ${FILE} | \
	  awk '{ temp = $1":"$2"[Amel4-5]"; $3 = temp; print $1"\t"$3"\t"$S2/28700"\t"$2 }' - \
	  >  ${DUMP}/selection/ehh/CHR/${CHR}.map

	# CHR cleanup
	rm ${FILE}
	rm ${DUMP}/selection/ehh/IDs

done



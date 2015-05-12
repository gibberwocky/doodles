#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# Set-up
module load bioinfo/bcftools

# ==============================================================================
# haps.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Haplotype alleles\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -v <VCF file> -o <output path> -p <population name> -s <VCF column start> -n <number of columns> [options]\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-x <T/F>\tInclude sites with missing genotypes [default = T]\e[0m"
  echo -e "\e[92m\t-h <T/F>\tInclude sites with heterozygous genotypes [default = F]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script creates hap and map files for a population within a VCF file.\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mBCFtools\e[0m"
  echo -e "\n"

}

# First 9 columns in VCF are annotation
# POP="JFM" START=10 N=30
# POP="OTH" START=40 N=30

VCF=
OUT=
POP=
START=
N=
MISSING=T
HET=F

while getopts ":v:o:p:s:n:x:h:" opt; do
  case $opt in
    v) VCF=${OPTARG};;
    o) OUT=${OPTARG};;
    p) POP=${OPTARG};;
    s) START=${OPTARG};;
    n) N=${OPTARG};;
    x) MISSING=${OPTGARG};;
    h) HET=${OPTARG};;
  esac
done

if [[ -z ${VCF} ]] | [[ -z ${OUT} ]] | [[ -z ${POP} ]] | [[ -z ${START}	]] | [[	-z ${N} ]]
then
  usage
  exit 1
fi



# Input VCF file
FILE=${OUT}/${POP}-tmp.vcf

#CHRS=('NC_007070.3' 'NC_007071.3' 'NC_007072.3' 'NC_007073.3' \
#  'NC_007074.3' 'NC_007075.3' 'NC_007076.3' 'NC_007077.3' 'NC_007078.3' \
#  'NC_007079.3' 'NC_007080.3' 'NC_007081.3' 'NC_007082.3' 'NC_007083.3' \
#  'NC_007084.3' 'NC_007085.3' 'NC_001566.1')

CHRS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16')


# Chromosomal loop
function join { local IFS="$1"; shift; echo "$*"; }

# Make directories for population and map files
mkdir -p ${OUT}/${POP}
mkdir -p ${OUT}/CHR

for CHR in ${CHRS[@]}
do

	# Extract chromosomes, exclude sites with uncalled genotype, drop header
	if [ "${MISSING}" = "T" ]; then
          if [ "${HET}" = "F" ]; then
 	    bcftools view -H -g ^het -O v --regions ${CHR} -o ${FILE} ${VCF}
	  fi
	fi
	if [ "${MISSING}" = "F" ]; then
       	  if [ "${HET}"	= "F" ]; then
	    bcftools view -H -g ^het -g ^miss -O v --regions ${CHR} -o ${FILE} ${VCF}
	  fi
	fi
        if [ "${MISSING}" = "T" ]; then
       	  if [ "${HET}"	= "T" ]; then
            bcftools view -H -O v --regions ${CHR} -o ${FILE} ${VCF}
	  fi
        fi
	if [ "${MISSING}" = "F" ]; then
       	  if [ "${HET}"	= "T" ]; then
            bcftools view -H -g ^miss -O v --regions ${CHR} -o ${FILE} ${VCF}
	  fi
        fi

	# Extract header and save sample IDs
	bcftools view -h  -O v -o ${OUT}/header-${POP}.vcf ${VCF}
#	tail --l 1 ${OUT}/header-${POP}.vcf > ${OUT}/IDs.${POP}
	rm ${OUT}/header-${POP}.vcf

	# Number of samples
	HAPFILE=${OUT}/${POP}/${POP}-${CHR}.hap
	rm -f ${HAPFILE}
	touch ${HAPFILE}

	for((I=$((${START}));I<=$((${START}+${N}-1));I++));
	do 
	  # Cut n paste
	  cut -f${I} ${FILE} > ${OUT}/test.${I}.${POP}
	  cut -c 1 ${OUT}/test.${I}.${POP} > ${OUT}/test.${I}a.${POP}
	  readarray tmp < ${OUT}/test.${I}a.${POP}
	  # Append data
	  echo ${tmp[@]} >> ${HAPFILE}
	  # Clean up
	  rm ${OUT}/test.${I}.${POP}
	  rm ${OUT}/test.${I}a.${POP}
	done


	# Extract mapfile details
	# 27027	= 1000000/37 (37cM/Mb)
	if [ -e "${OUT}/CHR/${CHR}.map" ]
	then
		echo "${OUT}/CHR/${CHR}.map exists, not re-writing"
	else
		cut -f 1-3 ${FILE} | awk '{ temp = $1":"$2"[Amel4-5]"; $3 = temp; print $1"\t"$3"\t"$2/27027"\t"$2 }' - >  ${OUT}/CHR/${CHR}.map
	fi

	# CHR cleanup
	rm ${FILE}
#	rm ${OUT}/IDs.${POP}

done



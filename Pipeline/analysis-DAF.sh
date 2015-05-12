#!/bin/bash
# IEKNGHER1ujm

# Set-up
module load bioinfo/bcftools
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"

# Root file name for Plink/Admixture input/output
DAF=${DUMP}/selection/daf
SITES=${DUMP}/vcfs/HARP-JFM-OTH



# ==============================================================================
# Repeat for each population's VCF file
# ==============================================================================

INFILE=OTH_beagle.vcf.gz
FILE=S2_JFM

# Reduce to chromosomes 1-16 + Mt
bcftools view --types snps -M 2 -O v \
  --regions 'NC_007070.3','NC_007071.3','NC_007072.3','NC_007073.3',\
'NC_007074.3','NC_007075.3','NC_007076.3','NC_007077.3','NC_007078.3',\
'NC_007079.3','NC_007080.3','NC_007081.3','NC_007082.3','NC_007083.3',\
'NC_007084.3','NC_007085.3','NC_001566.1' \
  -o ${DAF}/${FILE}.vcf \
  ${SITES}/${INFILE}
# Replace underscores in sample IDs
sed -i '/^#CHROM/s/_/-/g' ${DAF}/${FILE}.vcf
# Convert chromosome accessions to numbers
${PIPE}/accession-to-chr.sh -i ${DAF}/${FILE}.vcf


# 1) Delete header lines in VCF file, those begining with #
sed '/^#/ d' ${DAF}/${FILE}.vcf > ${DAF}/${FILE}.tmp1
# 2) Extract Chromosome and position
cut -f 1,2 ${DAF}/${FILE}.tmp1 > ${DAF}/${FILE}.tmp2
# 3) Extract AF from INFO field (method depends on format of INFO column in VCF file)
cut -f 8 ${DAF}/${FILE}.tmp1 > ${DAF}/${FILE}.tmp3
cut -d ';' -f 3 ${DAF}/${FILE}.tmp3 > ${DAF}/${FILE}.tmp4
#egrep -o 'AF=[^\t]+' ${DAF}/${FILE}.tmp1 > ${DAF}/${FILE}.tmp3
# 4) Remove "AF=" text
sed -i 's/AF=//g' ${DAF}/${FILE}.tmp4
# 5) Stick the annotation and AF values together
paste ${DAF}/${FILE}.tmp2 ${DAF}/${FILE}.tmp4 > ${DAF}/${FILE}.daf
# 6) Remove surplous files
rm ${DAF}/${FILE}.tmp*




# ==============================================================================
# Combine data to calculate DAF
# ==============================================================================

# Name of files
INFILE1=S2_JFM.daf
INFILE2=S2_OTH.daf
OUTFILE=S2_JFM-OTH.daf

# Join two datasets together
awk '{C[i=$1 OFS $2]; if(NR==FNR)A[i]=$3; else B[i]=$3} END{for(i in C) {print i,i in A?A[i]:".", i in B?B[i]:"."}}' OFS='\t' ${DAF}/${INFILE1} ${DAF}/${INFILE2} | sort > ${DAF}/${OUTFILE}
rm ${INFILE1}
rm ${INFILE2}

# Append column containing AFp1 - AFp2
awk '{ print $3-$4 }' ${DAF}/${OUTFILE} > ${DAF}/${OUTFILE}.tmp1
# Calculate mean
tmpm=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' ${DAF}/${OUTFILE}.tmp1)
# Calculate standard deviation
tmpd=$(awk '{sum+=$1; sumsq+=$1*$1} END {print sqrt(sumsq/NR - (sum/NR)**2)}' ${DAF}/${OUTFILE}.tmp1)
# Z-transformed tmp
awk -v VAR1=${tmpm} -v VAR2=${tmpd} '{ print ($1 - VAR1) / VAR2   }' ${DAF}/${OUTFILE}.tmp1 > ${DAF}/${OUTFILE}.tmp2
# Stick results together
paste ${DAF}/${OUTFILE} ${DAF}/${OUTFILE}.tmp1 ${DAF}/${OUTFILE}.tmp2 > ${DAF}/${OUTFILE}.tmp3

# Clean up
rm ${DAF}/${OUTFILE}
rm ${DAF}/${OUTFILE}.tmp1
rm ${DAF}/${OUTFILE}.tmp2
mv ${DAF}/${OUTFILE}.tmp3 ${DAF}/${OUTFILE}

# Output summary values to check max, min and average or "normal"
awk 'NR == 1 { max=$6; min=$6; sum=0 }
   { if ($6>max) max=$6; if ($6<min) min=$6; sum+=$6;}
   END {printf "Min: %d\tMax: %d\tAverage: %f\n", min, max, sum/NR}' ${DAF}/${OUTFILE}





# ==============================================================================
# Run DAF post-processing script in R to average results in windows
# ==============================================================================




#!/bin/bash
# IEKNGHER1ujm

# =============================================================================
# BreakDancer 
# gmt.genome.wustl.edu/packages/breakdancer/documentation.html
# =============================================================================
# Requires histogram module installing for Perl
# CPAN, install GD/Graph/histogram.pm

# BreakDancer's output file (.ctx) consists of the following columns:

# 1. Chromosome (BAM1)
# 2. Position (BAM1)
# 3. Orientation (BAM1)
# 4. Chromosome (BAM2)
# 5. Position (BAM2)
# 6. Orientation (BAM2)
# 7. Type of a SV
# 8. Size of a SV
# 9. Confidence Score
# 10. Total number of supporting read pairs
# 11. Total number of supporting read pairs from each map file
# 12. Estimated allele frequency (?)
# 13. Software version (?)
# 14. The run parameters (?)

# Columns 1-3 and 4-6 are used to specify the coordinates of the SV breakpoints
#	The orientation is a string that records the number of reads mapped to 
#	the plus (+) or the minus (-) strand in the anchoring regions.
# Column 7 is the type of SV detected
#	DEL (deletions), INS (insertion), INV (inversion), 
#	ITX (intra-chromosomal translocation), 
#	CTX (inter-chromosomal translocation), and Unknown. 
# Column 8 is the size of the SV in bp (meaningless for SV type ITX)
# Column 9 is the confidence score associated with the prediction. 
# Column 11 can be used to dissect the origin of the supporting read pairs
# 	For example, one may want to give SVs that are supported by more than 
#	one library higher confidence than those detected in only one library.  
#	It can also be used to distinguish somatic events from the germline, 
#	i.e., those detected in only the tumor libraries versus those detected 
#	in both the tumor and the normal libraries.
# Column 12 is currently a placeholder for displaying estimated allele frequency. 
#	The allele frequencies estimated in this version are not accurate and 
#	should not be trusted.
# Column 13 and 14 are information useful to reproduce the results.



# Pre-reqs
module load bioinfo/Java7
BDANCE="/usr/local/bioinfo/src/BreakDancer/current/bin"
BAM2CFG="/usr/local/bioinfo/src/BreakDancer/current/perl"
DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
AMEL="/save/seqapipop/Data/Apis-mellifera"

# Cretea file containing list of paths and file names of BAMs to merge
POP="OTH"
BAMS=(${AMEL}/SeqApiPop/${POP}/*/*_bootstrap.bam)
FILE=${DUMP}/breakdancer/${POP}.bams
for ((i=0;i<${#BAMS[@]};i++));
do
  if [[ $i = 0 ]]
    then printf "${BAMS[$i]}\n" > "${FILE}"
  fi
  if [[ $i > 0 ]]
    then printf "${BAMS[$i]}\n" >> "${FILE}"
  fi
done

# Merge BAMs in file using samtools, sort and index
POP="JFM"
qsub -q unlimitq -l mem=4G -l h_vmem=32G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/mergeBAMs-v2.sh \
      -n ${POP} \
      -b ${DUMP}/breakdancer/${POP}.bams \
      -o ${DUMP}/breakdancer

# Breakdancer outputs
OUT="JFM-OTH"

# Run BreakDancer
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/breakdancer.sh -a ${DUMP}/breakdancer/JFM.bam \
      -b ${DUMP}/breakdancer/OTH.bam -o ${DUMP}/breakdancer -s ${OUT}


# Extract data from resulting ctx file
# For instance, SVs with confidence score 99+
awk '$9 >=99' ${DUMP}/breakdancer/${OUT}.ctx \
  > ${DUMP}/breakdancer/${OUT}_conf99.ctx
# Change NCBI accessions to Ensembl
qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs \
    -e ${DUMP}/logs \
    ${PIPE}/accession-to-chr.sh -i ${DUMP}/breakdancer/${OUT}_conf99.ctx 
# Remove unplaced contigs
sed -i '/^NW_/d' ${DUMP}/breakdancer/${OUT}_conf99.ctx 
# Separate each of the SV types
awk '$7 ~/DEL/' ${DUMP}/breakdancer/${OUT}_conf99.ctx \
  > ${DUMP}/breakdancer/${OUT}_conf99-del.ctx
awk '$7 ~/INS/' ${DUMP}/breakdancer/${OUT}_conf99.ctx \
  > ${DUMP}/breakdancer/${OUT}_conf99-ins.ctx
awk '$7 ~/INV/' ${DUMP}/breakdancer/${OUT}_conf99.ctx \
  > ${DUMP}/breakdancer/${OUT}_conf99-inv.ctx
awk '$7 ~/ITX/' ${DUMP}/breakdancer/${OUT}_conf99.ctx \
  > ${DUMP}/breakdancer/${OUT}_conf99-itx.ctx
awk '$7 ~/CTX/' ${DUMP}/breakdancer/${OUT}_conf99.ctx \
  > ${DUMP}/breakdancer/${OUT}_conf99-ctx.ctx
# To extract any SVs within a particular interval
chr=10
pos1=387001
pos2=393000
awk -v VAR1=${chr} -v VAR2=${pos1} -v VAR3=${pos2} \
  '$1 ==VAR1 && $2 >=VAR2 && $2 <=VAR3' \
  ${DUMP}/breakdancer/${OUT}_conf99.ctx \
  > ${DUMP}/breakdancer/${chr}-${pos1}-${pos2}.ctx
# Fst hits checked, nothing special to report





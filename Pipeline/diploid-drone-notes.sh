#!/bin/bash
# IEKNGHER1ujm

# =============================================================================
# Create diploid drones
# =============================================================================

module load bioinfo/bcftools

DUMP="/home/dwragg/work/Analysis"
PIPE="/save/seqapipop/Scripts"
PERL=${PIPE}/Perl
#SITES=${DUMP}/vcfs/HARP-JFM-OTH-07112014
SITES=${DUMP}/vcfs/plink


# Identify samples for diploidy
POP="JFM"
SUFFIX="-33k_beagle"
readarray HAPS < ${DUMP}/vcfs/lists/${POP}_samples.list

# Create file to record sample IDs for merging
rm -f ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
touch ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
# Loop through all pairs
x=0
for ((i=0; i<${#HAPS[@]}; i=i+2))
do
  # Write sample ID to file
  printf "${SITES}/%s.vcf.gz\n" "${POP}-DIP${x}" \
    >> ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
  # Create directory for logs
  mkdir -p ${DUMP}/logs/${POP}-DIP${x}
  # Create diploid from samples i and i+1
  qsub -q unlimitq -l mem=4G -l h_vmem=16G -pe parallel_smp 1 \
    -o ${DUMP}/logs/${POP}-DIP${x} \
    -e ${DUMP}/logs/${POP}-DIP${x} \
    ${PIPE}/diploid-drone.sh \
      -v ${SITES}/${POP}${SUFFIX}.vcf.gz \
      -a ${HAPS[${i}]} \
      -b ${HAPS[${i}+1]} \
      -s ${POP}-DIP${x} \
      -o ${SITES}
  x=$((x+1))
done

# Merge diploids into single VCF
bcftools merge -m snps -O z \
  -l ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list \
  -o ${SITES}/${POP}_diploid.vcf.gz
tabix -p vcf ${SITES}/${POP}_diploid.vcf.gz
# Remove unnecessary files
readarray FILES < ${DUMP}/vcfs/lists/${POP}-diploid_vcfs.list
rm ${FILES[@]}
tmp=$(printf "%s.tbi " ${FILES[@]})
rm ${tmp}







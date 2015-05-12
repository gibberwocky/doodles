#!/bin/bash
# IEKNGHER1ujm

# M lineage samples
M_Harpur=('SRR957058' 'SRR957059' 'SRR957060' 'SRR957061' 'SRR957062' 'SRR957063' 'SRR957064' 'SRR957089' 'SRR957090')
M_Wallberg=('SRS549668' 'SRS549718' 'SRS549719' 'SRS549720' 'SRS549721' 'SRS549722' 'SRS549723' 'SRS549760' 'SRS549761' 'SRS549762' 'SRS549763')
M=(${M_Harpur[@]})
#M=(${M_Harpur[@]} ${M_Wallberg[@]})

# O lineage samples
O_Harpur=('SRR957065' 'SRR957066' 'SRR957067' 'SRR957068' 'SRR957069' 'SRR957070' 'SRR957071' 'SRR957072' 'SRR957073' 'SRR957074')
O_Wallberg=('SRS549663' 'SRS549664' 'SRS549665' 'SRS549666' 'SRS549667' 'SRS549713' 'SRS549714' 'SRS549715' 'SRS549716' 'SRS549717')
O=(${O_Harpur[@]})
#O=(${O_Harpur[@]} ${O_Wallberg[@]})

# A lineage samples
A_Harpur=('SRR957075' 'SRR957076' 'SRR957077' 'SRR957078' 'SRR957091' 'SRR957092' 'SRR957093' 'SRR957094' 'SRR957095' 'SRR957096' 'SRR957097')
A_Wallberg=('SRS549768' 'SRS549778')
A=(${A_Harpur[@]})
#A=(${A_Harpur[@]} ${A_Wallberg[@]})

# C lineage samples
C_Harpur=('SRR957080' 'SRR957081' 'SRR957082' 'SRR957083' 'SRR957084' 'SRR957085' 'SRR957086' 'SRR957087' 'SRR957088')
C_Wallberg=('SRS549629' 'SRS549630' 'SRS549631' 'SRS549632' 'SRS549635' 'SRS549685' 'SRS549686' 'SRS549687')
C_Liu=('SRR1425485' 'SRR1425486' 'SRR1425470' 'SRR1425471' 'SRR1425450' 'SRR1425453')
C=(${C_Harpur[@]})
#C=(${C_Harpur[@]} ${C_Wallberg[@]})


# Location of soft-linked VCF files
AMEL=/work/dwragg/Apis-mellifera

# M lineage VCF and BAM arrays
tmp=(${M[@]/#/"${AMEL}/*/vcfs/"})
VCF_M=(${tmp[@]/%/'_clean.vcf.gz'})
tmp=(${M[@]/#/"${AMEL}/*/bams/"})
BAM_M=(${tmp[@]/%/'_bootstrap.bam'})

# A lineage VCF and BAM arrays
tmp=(${A[@]/#/"${AMEL}/*/vcfs/"})
VCF_A=(${tmp[@]/%/'_clean.vcf.gz'})
tmp=(${A[@]/#/"${AMEL}/*/bams/"})
BAM_A=(${tmp[@]/%/'_bootstrap.bam'})

# C lineage VCF and BAM arrays
tmp=(${C[@]/#/"${AMEL}/*/vcfs/"})
VCF_C=(${tmp[@]/%/'_clean.vcf.gz'})
tmp=(${C[@]/#/"${AMEL}/*/bams/"})
BAM_C=(${tmp[@]/%/'_bootstrap.bam'})

# O lineage VCF and BAM arrays
tmp=(${O[@]/#/"${AMEL}/*/vcfs/"})
VCF_O=(${tmp[@]/%/'_clean.vcf.gz'})
tmp=(${O[@]/#/"${AMEL}/*/bams/"})
BAM_O=(${tmp[@]/%/'_bootstrap.bam'})

# Print list of Harpur samples
printf "%s\n" "${A_Harpur[@]}" > ${DUMP}/vcfs/lists/A_Harpur_samples.list
printf "%s\n" "${M_Harpur[@]}" > ${DUMP}/vcfs/lists/M_Harpur_samples.list
printf "%s\n" "${C_Harpur[@]}" > ${DUMP}/vcfs/lists/C_Harpur_samples.list
printf "%s\n" "${O_Harpur[@]}" > ${DUMP}/vcfs/lists/O_Harpur_samples.list

# Print list of Harpur + Wallberg samples
printf "%s\n" "${A[@]}" > ${DUMP}/vcfs/lists/A_HW_samples.list
printf "%s\n" "${M[@]}" > ${DUMP}/vcfs/lists/M_HW_samples.list
printf "%s\n" "${C[@]}" > ${DUMP}/vcfs/lists/C_HW_samples.list
printf "%s\n" "${O[@]}" > ${DUMP}/vcfs/lists/O_HW_samples.list



POP="BR"
cd /save/seqapipop/Data/Apis-mellifera/SeqApiPop/${POP}
SAMPLES=(*)
printf "%s\n" "${SAMPLES[@]}" > ${DUMP}/vcfs/lists/${POP}_samples.list









# ==========================================================
# Download SRA data and convert to FASTQ
# ==========================================================
PIPE="/save/seqapipop/Scripts"
DATA=('SRR1425486' 'SRR1425470' 'SRR1425471' 'SRR1425450' 'SRR1425453')

OUT=/work/dwragg/SRA
for SAMPLE in ${DATA[@]}
do
  mkdir -p ${OUT}/${SAMPLE}/logs
  qsub -q workq -l mem=8G -l h_vmem=8G -pe parallel_smp 8 \
    -o ${OUT}/${SAMPLE}/logs \
    -e ${OUT}/${SAMPLE}/logs \
     ${PIPE}/SRA-fastq.sh \
      -s ${SAMPLE} \
      -o ${OUT}
done

DATA=('SRR1425486' 'SRR1425470' 'SRR1425471' 'SRR1425450' 'SRR1425453')
for SAMPLE in ${DATA[@]}
do
  rm ${OUT}/${SAMPLE}/${SAMPLE}_R*.fastq.gz
done

# Requires $SAMPLE_R1.fastq.gz $SAMPLE_R2.fastq.gz in $OUT/$SAMPLE folder
EMBOSS="/usr/local/bioinfo/src/EMBOSS/EMBOSS-6.4.0/bin"
DATA=('SRR1425486' 'SRR1425470' 'SRR1425471' 'SRR1425450' 'SRR1425453')
for SAMPLE in ${DATA[@]}
do
#  cp ${OUT}/${SAMPLE}/${SAMPLE}_1.fastq.gz ${OUT}/${SAMPLE}/${SAMPLE}_R1.fastq.gz
#  cp ${OUT}/${SAMPLE}/${SAMPLE}_2.fastq.gz ${OUT}/${SAMPLE}/${SAMPLE}_R2.fastq.gz
  qsub -q workq -l mem=8G -l h_vmem=8G -pe parallel_smp 8 \
    -o ${OUT}/${SAMPLE}/logs \
    -e ${OUT}/${SAMPLE}/logs \
  ${PIPE}/phred64-2-phred22.sh -e ${EMBOSS} -s ${SAMPLE} -f ${OUT} -p "T"
done

# Above step needs running once fastq files available



# C_LIU samples have extremely high quality scores, SRR1425453 to be re-done
# map.sh	stat.sh		snps.sh
DATA=('SRR1425485' 'SRR1425486' 'SRR1425470' 'SRR1425471' 'SRR1425450' 'SRR1425453')
for SAMPLE in ${DATA[@]}
do
  qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${OUT}/${SAMPLE}/logs \
    -e ${OUT}/${SAMPLE}/logs \
    ${PIPE}/snps.sh -s ${SAMPLE} -f ${OUT}/${SAMPLE} -o ${OUT}
done






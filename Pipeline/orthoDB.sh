#!bin/bash
# IEKNGHER1ujm

# BASH
# 	1) Amel_RNA BLASTX Dros_PEP
# 	2) Dros_RNA BLASTX Amel_PEP
# R
#	1) Filter out alignments < 50 in length
#	2) Sort by transcript ID and decreasing e-value
#	3) Remove duplicates 
#		(in other words, steps 2 and 3 result in the alignment with the
#		lowest e-value being retained)
# 	4) Cross-reference RNA and transcript IDs with gene IDs for each species
#	5) Merge the two tables on matching Gene IDs for each species

# Result 
#	16713 orthologous transcripts representing 6927 genes



# ApisDros orthoDB
DUMP="/work/dwragg/Analysis/orthoDB"
PIPE="/home/dwragg/work/Pipeline"

# Ensembl transcript data
TRANSCRIPTS="/save/dwragg/Transcripts"
Amel_RNA="Apis_mellifera.GCA_000002195.1.23.cdna.all.fa.gz"
Amel_PEP="Apis_mellifera.GCA_000002195.1.23.pep.all.fa.gz"
Dros_RNA="Drosophila_melanogaster.BDGP5.23.cdna.all.fa.gz"
Dros_PEP="Drosophila_melanogaster.BDGP5.23.pep.all.fa.gz"

# Extract to work partition
gunzip -c ${TRANSCRIPTS}/${Amel_RNA} > ${TRANSCRIPTS}//Amel_RNA.fa
gunzip -c ${TRANSCRIPTS}/${Amel_PEP} > ${TRANSCRIPTS}//Amel_PEP.fa
gunzip -c ${TRANSCRIPTS}/${Dros_RNA} > ${TRANSCRIPTS}//Dros_RNA.fa
gunzip -c ${TRANSCRIPTS}/${Dros_PEP} > ${TRANSCRIPTS}//Dros_PEP.fa

# Format BLAST databases
formatdb -i ${TRANSCRIPTS}/Amel_RNA.fa -o T -p F
formatdb -i ${TRANSCRIPTS}/Dros_RNA.fa -o T -p F
formatdb -i ${TRANSCRIPTS}/Amel_PEP.fa -o T -p T
formatdb -i ${TRANSCRIPTS}/Dros_PEP.fa -o T -p T

# Count number of transcripts (27k)
#grep -e "^>" Amel_RNA.fa | wc -l


# ==============================================================================
# For each RNA set BLASTX to opposing PEP set
# ==============================================================================
# -i input fasta (eg Amel_RNA)
# -d database fasta (eg Dros_PEP)
# -t sets the program to use
# -f sets the output format (9 = tabular with comments)
# -e sets the e-value threshold (Sweden used 0.5)
# -m sets to use megablast
# -o output file for stdout

qsub -q unlimitq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
  -o ${DUMP} \
  -e ${DUMP} \
  ${PIPE}/orthologs.sh -i ${TRANSCRIPTS}/Dros_RNA.fa -d ${TRANSCRIPTS}/Amel_PEP.fa \
    -t blastx -f 8 -e 0.5 -m T -o ${DUMP}/Dros_RNA-Amel_PEP.txt

# ==============================================================================
# Add field names to start of BLASTX results and remove comment lines
# ==============================================================================
printf 'Query_ID\tSubject_ID\tpc_ID\tALIGN_LEN\tn_MISMATCH\tn_GAPs\tQ_Start\tQ_End\tS_Start\tS_End\te-Val\tBit_Score\n' > ${DUMP}/fields.txt
cat ${DUMP}/fields.txt ${DUMP}/Amel_RNA-Dros_PEP.txt > Amel_RNA-Dros_PEP
cat ${DUMP}/fields.txt ${DUMP}/Dros_RNA-Amel_PEP.txt > Dros_RNA-Amel_PEP
sed -i '/^#/ d' Amel_RNA-Dros_PEP
sed -i '/^#/ d' Dros_RNA-Amel_PEP
sed -i 's/lcl|//g' Amel_RNA-Dros_PEP
sed -i 's/lcl|//g' Dros_RNA-Amel_PEP
rm *PEP.txt

# The gb accessions belong to the GENPEPT database (DAVID)
# To use in gProfiler first need to use gOrth to get mouse ortholog ensembl IDs



# To merge the two tables, the column names must match in addition to some
# 	of the column values. This may be a problem when looking at different
#	types of transcript (e.g. RNA vs peptide)
#	Will need to cross-reference transcripts IDs I think...
# The following extracts the transcript and gene IDs from the fasta files so that
# they can be used as an index for cross-referencing transcript IDs
ID="Amel_PEP"
#ID="Dros_PEP"
grep -e "^>" ${TRANSCRIPTS}/${ID}.fa > ${DUMP}/${ID}.tmp
sed -i 's/^>//g' ${DUMP}/${ID}.tmp
sed -i 's/gene://g' ${DUMP}/${ID}.tmp
sed -i 's/transcript://g' ${DUMP}/${ID}.tmp
cut -f 1,4,5 -d ' ' ${DUMP}/${ID}.tmp >  ${DUMP}/${ID}.db

# Merge Subject_ID from Amel_RNA-Dros_PEP (FBppxxxx) with Dros_PEP.db store FBtrxxx

# ==============================================================================
# Process output tables in R to get best alignment matches
# ==============================================================================

setwd("/work/dwragg/Analysis/orthoDB")

# Load peptide-transcript indexes
amel_pep_idx <- read.table("Amel_PEP.db", header=F, stringsAsFactors=F, row.names=NULL)
names(amel_pep_idx)[1:3] <- c("Subject_ID", "Gene", "Amel_Trans_ID")
dros_pep_idx <- read.table("Dros_PEP.db", header=F, stringsAsFactors=F, row.names=NULL)
names(dros_pep_idx)[1:3] <- c("Subject_ID", "Gene", "Dros_Trans_ID")


# load DROSOPHILA RNA -> APIS PEP blastx results
dat <- read.table("Dros_RNA-Amel_PEP", header=T, stringsAsFactors=F, row.names=NULL)
DROS <- merge(dat, amel_pep_idx, by="Subject_ID")
names(DROS)[2] <- "Dros_Trans_ID"
# Retain only alignments spanning min of 50 peptides
DROS <- DROS[which(DROS$ALIGN_LEN >=50),]
# Sort the table by Amel_Trans_ID (RNA transcript) and decreasing e-value
tmp <- DROS[order(DROS$Dros_Trans_ID, DROS$e.Val, decreasing=F ), ]
# Remove duplicate Amel_Trans_ID (in otherwords, retain that with lowest e-value)
DROS <- tmp[!duplicated(tmp$Dros_Trans_ID), ]
# Update with gene and peptide ID data
names(dros_pep_idx)[1:3] <- c("Dros_Pep_ID", "Dros_Gene_ID", "Dros_Trans_ID")
DROS <- merge(DROS, dros_pep_idx, by="Dros_Trans_ID")
names(DROS)[2] <- "Amel_Pep_ID"
names(DROS)[13] <- "Amel_Gene_ID"

# load APIS RNA -> DROSOPHILA PEP blastx results
names(dros_pep_idx)[1:3] <- c("Subject_ID", "Gene", "Dros_Trans_ID")
dat <- read.table("Amel_RNA-Dros_PEP", header=T, stringsAsFactors=F, row.names=NULL)
AMEL <- merge(dat, dros_pep_idx, by="Subject_ID")
names(AMEL)[2] <- "Amel_Trans_ID"
# Retain only alignments spanning min of 50 peptides
AMEL <- AMEL[which(AMEL$ALIGN_LEN >=50),]
# Sort the table by Amel_Trans_ID (RNA transcript) and decreasing e-value
tmp <- AMEL[order(AMEL$Amel_Trans_ID, AMEL$e.Val, decreasing=F ), ]
# Remove duplicate Amel_Trans_ID (in otherwords, retain that with lowest e-value)
AMEL <- tmp[!duplicated(tmp$Amel_Trans_ID), ]
# Update with gene and peptide ID data
names(amel_pep_idx)[1:3] <- c("Amel_Pep_ID", "Amel_Gene_ID", "Amel_Trans_ID")
AMEL <- merge(AMEL, amel_pep_idx, by="Amel_Trans_ID")
names(AMEL)[2] <- "Dros_Pep_ID"
names(AMEL)[13] <- "Dros_Gene_ID"


# <================= to here is good.


# Combine the two sets of results on matching Dros and Amel Gene IDs
dat <- merge(AMEL, DROS, by=c("Dros_Gene_ID", "Amel_Gene_ID"), suffixes=c(".AMEL", ".DROS"))

write.table(dat, "orthologs.txt", append=F, quote=F, sep="\t", row.names=F, col.names=T)
write.table(AMEL, "AMEL.txt", append=F,  quote=F, sep="\t", row.names=F, col.names=T)
write.table(DROS, "DROS.txt", append=F,  quote=F, sep="\t", row.names=F, col.names=T)



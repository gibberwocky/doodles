#!/usr/bin/env Rscript

args <- commandArgs(T)

# (1) sample
# (2) output
# (3) platypus_pass.vcf.gz
# (4) GATK_UG_pass.vcf.gz
# (5) pileup_hom.vcf.gz
# (6) BAYSIC.vcf.gz


# VennDiagram library
library("VennDiagram")

# Read GATK file
tmp <- read.table(args[4], sep="\t", stringsAsFactors=F)
# Combine chromosome and position elements into a vector
varA <- paste(tmp[,1], tmp[,2], sep="-")
# Create a vector of equal length containing category name
nomA <- rep("GATK", length(varA))
rm(tmp)

# Read Mpileup file
tmp <- read.table(args[5], sep="\t", stringsAsFactors=F)
# Combine chromosome and position elements into a vector
varB <- paste(tmp[,1], tmp[,2], sep="-")
# Create a vector of equal length containing category name
nomB <- rep("Pileup", length(varB))
rm(tmp)

# Read Platypus file
tmp <- read.table(args[3], sep="\t", stringsAsFactors=F)
# Combine chromosome and position elements into a vector
varC <- paste(tmp[,1], tmp[,2], sep="-")
# Create a vector of equal length containing category name
nomC <- rep("Platypus", length(varC))
rm(tmp)

# Read BAYSIC file
tmp <- read.table(args[6], sep="\t", stringsAsFactors=F)
# Combine chromosome and position elements into a vector
varD <- paste(tmp[,1], tmp[,2], sep="-")
# Create a vector of equal length containing category name
nomD <- rep("BAYSIC", length(varD))
rm(tmp)


# Plot Venn diagram using VennDiagram package
venn.diagram(list(GATK=varA, Pileup=varB, Platypus=varC, BAYSIC=varD),
	fill=c("red", "green", "blue", "yellow"),
	alpha=c(0.5, 0.5, 0.5, 0.5), cex=0.8,
	cat.fontfamily="sans", cat.fontface=2,
	fontfamily="sans", fontface=1,
	lty=1, euler.d=TRUE, scaled=TRUE,
	reverse=TRUE, compression='lzw',
	resolution=300, height=5, width=5, units='in', 
	filename=args[2]);



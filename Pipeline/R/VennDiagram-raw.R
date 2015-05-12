#!/usr/bin/env Rscript

args <- commandArgs(T)

# VennDiagram library
library("VennDiagram")

# Read A
tmp <- read.table(args[3], sep="\t", stringsAsFactors=F)
# Combine chromosome and position elements into a vector
varA <- paste(tmp[,1], tmp[,2], sep="-")
# Create a vector of equal length containing category name
nomA <- rep("A", length(varA))
rm(tmp)

# Read M
tmp <- read.table(args[4], sep="\t", stringsAsFactors=F)
# Combine chromosome and position elements into a vector
varB <- paste(tmp[,1], tmp[,2], sep="-")
# Create a vector of equal length containing category name
nomB <- rep("M", length(varB))
rm(tmp)

# Read C
tmp <- read.table(args[5], sep="\t", stringsAsFactors=F)
# Combine chromosome and position elements into a vector
varC <- paste(tmp[,1], tmp[,2], sep="-")
# Create a vector of equal length containing category name
nomC <- rep("C", length(varC))
rm(tmp)

# Read O
tmp <- read.table(args[6], sep="\t", stringsAsFactors=F)
# Combine chromosome and position elements into a vector
varD <- paste(tmp[,1], tmp[,2], sep="-")
# Create a vector of equal length containing category name
nomD <- rep("O", length(varD))
rm(tmp)

# Read JFM-OTH raw
tmp <- read.table(args[7], sep="\t", stringsAsFactors=F)
# Combine chromosome and position elements into a vector
varE <- paste(tmp[,1], tmp[,2], sep="-")
# Create a vector of equal length containing category name
nomE <- rep("Haploid", length(varE))
rm(tmp)



# Plot Venn diagram using VennDiagram package
venn.diagram(list(A=varA, M=varB, C=varC, O=varD, Haploid=varE),
	fill=c("red", "green", "blue", "yellow", "purple"),
	alpha=c(0.5, 0.5, 0.5, 0.5, 0.5), cex=0.5, cat.cex=0.6,
	cat.fontfamily="sans", cat.fontface=2,
	fontfamily="sans", fontface=1,
	lty=1, euler.d=TRUE, scaled=TRUE,
	reverse=TRUE, compression='lzw', margin=0.1,
	resolution=300, height=6, width=6, units='in', 
	filename=args[2]);



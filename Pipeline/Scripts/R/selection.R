#!/usr/bin/env Rscript

# ==============================================================================
# NOTES
# ==============================================================================

# The following functions are inlcuded in this script:

# func_majA	- To identify and count major allele
# func_minA	- To identify and count minor allele

# ==============================================================================

args <- commandArgs(T)
# [1] path
# [2] ped
# [3] map
# [4] pheno
# [5] chromosome
# [6] window size
# [7] chrinfo
# [8] output

#args <- c("/home/dwragg/work/Analysis/vcfs", "EurBee-bayes.ped", "EurBee-R.map", "EurBee-bayes.pheno", "1", "25000", "/save/dwragg/Apis/chrInfo", "/home/dwragg/work/Analysis/selection")

library("GenABEL")

# Set working directory
setwd(args[1])

# Import PED data
genFile <- paste(args[8], "/gen", args[5], sep="")
convert.snp.ped(pedfile=args[2], mapfile=args[3], out=genFile, format = "premakeped", traits = 1, strand = "u", bcast = 10000000, wslash=F, mapHasHeaderLine=F)
raw <- load.gwaa.data(phe=args[4], gen = genFile, force = T, sort=T)

# Load table of chromosome lengths
chrInfo <- read.table(args[7], sep="\t", header=F, strings=F)



# ==============================================================================
# Generate windows
# ==============================================================================

# Get length
var_chr <- as.numeric(args[5])
var_chrL <- as.numeric(chrInfo[which(chrInfo$V1==var_chr), 2])

# Specify window size
var_winL <- as.numeric(args[6])
var_overlap <- 0.5

# Calculates intervals based on parameters
start <- as.integer(seq(1,var_chrL, by=var_winL*(var_overlap)))
chr <- rep(var_chr, length(start))
stop <- as.integer(start + var_winL)
bed <- NULL
bed <- data.frame(chr, start, stop, stringsAsFactors=F)
bed[nrow(bed),ncol(bed)] <- var_chrL



# ==============================================================================
# Count genotype and major/minor allele frequencies
# ==============================================================================

# Subset to chromosome
data <- raw[,which(raw@gtdata@chromosome==var_chr)]
rm(raw)

# Get list of SNP positions
var_Map <- data@gtdata@map

# Get list of SNP names
var_SNP <- data@gtdata@snpnames

# Get number of SNPs per window
nSNPs <- lapply(as.integer(rownames(bed)), 
	function(x) sum(var_Map%in%bed[x,]$start:bed[x,]$stop))

# Get list of SNPs per window
SNPs <- lapply(as.integer(rownames(bed)),
	function(x) var_SNP[var_Map%in%bed[x,]$start:bed[x,]$stop])

# Extract number per genotype
summarySNPs <- summary(data)
P11 <- lapply(SNPs, function(x) if(length(x)>0) {
	sum(summarySNPs[which(row.names(summarySNPs)%in%x),]$P.11)} else {0})
P12 <- lapply(SNPs, function(x) if(length(x)>0) {
	sum(summarySNPs[which(row.names(summarySNPs)%in%x),]$P.12)} else {0})
P22 <- lapply(SNPs, function(x) if(length(x)>0) {
	sum(summarySNPs[which(row.names(summarySNPs)%in%x),]$P.22)} else {0})

# Function for calculating major allele
func_majA <- function(x){
	tmp <- summary(data[,x])
	2 * (sum(unlist(lapply(rownames(tmp), function(y) max(tmp[y,c(9,11)]))))) +
	sum(unlist(lapply(rownames(tmp), function(y) tmp[y,10])))
}
# Function for calculating minor allele
func_minA <- function(x){
	tmp <- summary(data[,x])
	2 * sum(unlist(lapply(rownames(tmp), function(y) min(tmp[y,c(9,11)])))) +
	sum(unlist(lapply(rownames(tmp), function(y) tmp[y,10])))
}
# Run major allele function
aMj <- lapply(SNPs, function(x) if(length(x)>0) {func_majA(x)} else {0})
# Run minor allele function
aMn <- lapply(SNPs, function(x) if(length(x)>0) {func_minA(x)} else {0})


# Record values and label columns
bed[,4] <- unlist(nSNPs)
names(bed)[4] <- "n_SNPs"

bed[,5] <- unlist(P11)
names(bed)[5] <- "n_P11"

bed[,6] <- unlist(P12)
names(bed)[6] <- "n_P12"

bed[,7] <- unlist(P22)
names(bed)[7] <- "n_P22"

bed[,8] <- unlist(aMj)
names(bed)[8] <- "n_Maj"

bed[,9] <- unlist(aMn)
names(bed)[9] <- "n_Min"



# ==============================================================================
# Calculate pooled heterozygosity from major/minor allele frequencies
# ==============================================================================

# Calculate and record pooled heterozgosity = (2*aMj*aMn) / (aMj+aMn)^2
bed[,10] <- round(( 2 * bed[,8] * bed[,9] ) / ( bed[,8] + bed[,9] ) ^2, 3)
names(bed)[10] <- "Hp"

# Standardize Hp using mean (if normally distributed) or median (if skewed)
# Negative is sweep (positive selection)
# Positive is hypervariable (balancing selection)
var_Mean <- mean(bed$Hp, na.rm=T)
var_Median <- median(bed$Hp, na.rm=T)
var_Sdev <- sd(bed$Hp, na.rm=T)
bed[,11] <- round((bed[,10]-var_Mean)/var_Sdev, 3)
names(bed)[11] <- "Z_Hp"


# Write results to file
fname <- paste(args[8], "/Hp-chr", var_chr, "-win-", var_winL, ".txt", sep="")
write.table(bed, fname, sep="\t", row.names=T, col.names=T, append=F, quote=F)








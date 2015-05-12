#!/usr/bin/env Rscript

# ==============================================================================
# NOTES
# ==============================================================================

# The following functions are inlcuded in this script:

# func_majA	- To identify and count major allele
# func_minA	- To identify and count minor allele
# func_Fst 	- To calculate unbias Fst and Jost D using mean and median

# ==============================================================================

args <- commandArgs(T)
# [1] path
# [2] ped
# [3] map
# [4] pheno
# [5] chromosome
# [6] window size


library("GenABEL")

# Set working directory
setwd(args[1])

# Import PED data
convert.snp.ped(pedfile=args[2], mapfile=args[3], out="gen0i1.raw", format = "premakeped", traits = 1, strand = "u", bcast = 10000000, wslash=F, mapHasHeaderLine=F)
raw <- load.gwaa.data(phe=args[4], gen = "gen0i1.raw", force = T, sort=T)

# Load table of chromosome lengths
chrInfo <- read.table("EurBee.chrInfo", sep="\t", header=F, strings=F)

# ==============================================================================
# Generate windows
# ==============================================================================

# Get chromosome 24 Length
var_chr <- args[5]
var_chrL <- chrInfo[which(chrInfo$V1==var_chr), 2]

# Specify window size
var_winL <- args[6]
var_overlap <- 0.5

# Calculates intervals based on parameters
start <- as.integer(seq(1,var_chrL, by=var_winL*(var_overlap)))
chr <- rep(var_chr, length(start))
stop <- as.integer(start + var_winL)
bed <- NULL
bed <- data.frame(chr, start, stop, stringsAsFactors=F)
bed[nrow(bed),ncol(bed)] <- var_chrL

# ==============================================================================



# ==============================================================================
# Count genotype and major/minor allele frequencies
# ==============================================================================

# Subset to chromosome
data <- raw[,which(data@gtdata@chromosome==var_chr)]
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



# ==============================================================================
# Calculate pooled heterozygosity from major/minor allele frequencies
# ==============================================================================

# Calculate and record pooled heterozgosity = (2*aMj*aMn) / (aMj+aMn)^2
bed[,10] <- ( 2 * bed[,8] * bed[,9] ) / ( bed[,8] + bed[,9] ) ^2
names(bed)[10] <- "Hp"

# Check if normally distributed
pdf("Hp-dist.pdf", height=5, width=5)
hist(bed$Hp)
dev.off()

# Standardize Hp using mean (if normally distributed) or median (if skewed)
# Negative is sweep (positive selection)
# Positive is hypervariable (balancing selection)
var_Mean <- mean(bed$Hp, na.rm=T)
var_Median <- median(bed$Hp, na.rm=T)
var_Sdev <- sd(bed$Hp, na.rm=T)
bed[,11] <- round((bed[,10]-var_Mean)/var_Sdev, 2)
names(bed)[11] <- "Z_Hp"

# Plot Hp
zscore_x_max = max(bed$Z_Hp, na.rm=T) + 1
zscore_x_min = min(bed$Z_Hp, na.rm=T) - 1
pos_mid <- ((bed$start + bed$stop) / 2)
max_y <- max(bed$Z_Hp, na.rm=T)
max_x <- max(pos_mid)
xfiglim <- round(max_x/40000)
pdf(file="Hp-Plot.pdf", width=5, height=5)
profile<-plot(bed$Z_Hp, x=pos_mid, ylim=range(zscore_x_min:zscore_x_max), 
	pch=".", col="blue", type="l", axes=FALSE, ylab="Zscore", xlab="", 
	cex.lab=0.6)
title(main = "Plot", cex.main=0.7)
axis(2, line=0, at=seq(zscore_x_min, zscore_x_max, by=1), lwd=1, lwd.ticks=1, 
	las=2, cex.axis=0.6)
axis(1, line=0, at=seq(0, max_x, by=1000000), lwd=2, lwd.ticks=1, las=2, 
	col="white", col.ticks="black", cex.axis=0.6 )
axis(1, line=0, at=c(0, max_x), labels=c("",""), lwd.ticks=0, lwd=1, cex=0.4)
abline(h=c(zscore_x_min:zscore_x_max),v=NULL, lty=3)
# Polygon is a bit weird, should investigate this to see why
polygon(c(0,pos_mid), c(0,bed$Z_Hp), col="blue")
dev.off()

# ==============================================================================



# ==============================================================================
# Calculate Fst and Jost's D
# ==============================================================================

# Function: unbias.Fst
func_Fst <- function(dat.geno, pop1.inds, pop2.inds, SNPs){
	
    # Total population
    Ht.dat <- dat.geno[c(pop1.inds, pop2.inds),SNPs]
    Ht.11 <- summary(Ht.dat)$P.11
    Ht.12 <- summary(Ht.dat)$P.12
    Ht.22 <- summary(Ht.dat)$P.22
    Ht.n <- Ht.11 + Ht.12 + Ht.22
    Ht <- (((2*Ht.11)+Ht.12)/(2*Ht.n)) * (1 - (((2*Ht.11)+Ht.12)/(2*Ht.n))) * 2
    
    # Pop1
    Hs1.dat <- dat.geno[pop1.inds, SNPs]
    Hs1.11 <- summary(Hs1.dat)$P.11
    Hs1.12 <- summary(Hs1.dat)$P.12
    Hs1.22 <- summary(Hs1.dat)$P.22
    Hs1.n <- Hs1.11 + Hs1.12 + Hs1.22
    Hs1 <- (((2*Hs1.11)+Hs1.12)/(2*Hs1.n)) * (1 - (((2*Hs1.11)+Hs1.12)/(2*Hs1.n))) * 2
    
    # Pop2
    Hs2.dat <- dat.geno[pop2.inds, SNPs]
    Hs2.11 <- summary(Hs2.dat)$P.11
    Hs2.12 <- summary(Hs2.dat)$P.12
    Hs2.22 <- summary(Hs2.dat)$P.22
    Hs2.n <- Hs2.11 + Hs2.12 + Hs2.22
    Hs2 <- (((2*Hs2.11)+Hs2.12)/(2*Hs2.n)) * (1 - (((2*Hs2.11)+Hs2.12)/(2*Hs2.n))) * 2
    
    # Weighted Fst
    wHs <- ( (Hs1*Hs1.n) + (Hs2*Hs2.n) ) / Ht.n
    fst1 <- median((Ht-wHs)/Ht, na.rm=T)
    fst2 <- mean((Ht-wHs)/Ht, na.rm=T)

    # D (Joost, 2008)
    k <- 2
    D1 <- (k/(k-1)) * median ( (Ht-wHs) / (1-wHs), na.rm=T )
    D2 <- (k/(k-1)) * mean ( (Ht-wHs) / (1-wHs), na.rm=T )

    rm(k, dat.geno, pop1.inds, pop2.inds, SNPs, wHs)
    rm(Ht.dat, Ht.11, Ht.12, Ht.22, Ht.n, Ht)
    rm(Hs1.dat, Hs1.11, Hs1.12, Hs1.22, Hs1.n, Hs1)
    rm(Hs2.dat, Hs2.11, Hs2.12, Hs2.22, Hs2.n, Hs2)
    return(c(fst1, fst2, D1, D2))
  
}# End of function


# Apply function to data
fst <- lapply(SNPs, function(x) if(length(x)>0) {
	func_Fst(data@gtdata, 
		data@phdata[which(data@phdata$pop==1),]$id,
		data@phdata[which(data@phdata$pop==2),]$id,
		x)} else {c(0,0,0,0)})


bed[,12:15] <- unlist(fst)
names(bed)[12:15] <- c("Fst_median", "Fst_mean", "D_median", "D_mean")

# ==============================================================================


fname <- paste("selection-chr", var_chr, "-win-", var_winL, ".txt", sep="")
write.table(bed, fname, sep="\t", row.names=T, col.names=T, append=F, quote=F)








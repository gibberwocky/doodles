#!/usr/bin/env Rscript

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
# [9] overlap

#args <- c("/home/dwragg/work/Analysis/vcfs", "EurBee-bayes.ped", "EurBee-R.map", "EurBee-bayes.pheno", "1", "25000", "/save/dwragg/Apis/chrInfo", "/home/dwragg/work/Analysis/selection", "0.25")

library("GenABEL")
library("plyr")

# Set working directory
setwd(args[1])

# Import PED data
genFile <- paste(args[8], "/gen", args[5], sep="")
convert.snp.ped(pedfile=args[2], mapfile=args[3], out=genFile, format = "premakeped", traits = 1, strand = "u", bcast = 10000000, wslash=F, mapHasHeaderLine=F)
raw <- load.gwaa.data(phe=args[4], gen = genFile, force = T, sort=T)
qc0 <- check.marker(raw, maf = 0.01, odds=NA, p.lev = 0, ibs.mrk = 0 , het.fdr = 0, perid.call = 0, callrate = 0)
raw <- raw[qc0$idok, qc0$snpok]

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
var_overlap <- as.numeric(args[9])

# Calculates intervals based on parameters
#start <- as.integer(seq(1,var_chrL, by=var_winL*(var_overlap)))
start <- as.integer(seq(1,var_chrL, by=var_overlap))
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
cat("Get SNP positions\n")
var_Map <- data@gtdata@map

# Get list of SNP names
cat("Get SNP names\n")
var_SNP <- data@gtdata@snpnames

# Get number of SNPs per window
cat("Get number of SNPs per window\n")
#nSNPs <- lapply(as.integer(rownames(bed)), 
#	function(x) sum(var_Map%in%bed[x,]$start:bed[x,]$stop))
nSNPs <- llply(as.integer(rownames(bed)),
	function(x) sum(var_Map%in%bed[x,]$start:bed[x,]$stop),
	.progress="text", .inform=TRUE, .parallel=FALSE, .paropts=NULL)

# Get list of SNPs per window
cat("Get SNPs per window\n")
#SNPs <- lapply(as.integer(rownames(bed)),
#	function(x) var_SNP[var_Map%in%bed[x,]$start:bed[x,]$stop])
SNPs <- llply(as.integer(rownames(bed)),
       function(x) var_SNP[var_Map%in%bed[x,]$start:bed[x,]$stop],
	.progress="text", .inform=TRUE,	.parallel=FALSE, .paropts=NULL)


# Record values and label columns
bed[,4] <- unlist(nSNPs)
names(bed)[4] <- "n_SNPs"



# ==============================================================================
# Calculate Fst
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
    fst2 <- round(mean((Ht-wHs)/Ht, na.rm=T), 3)

    rm(dat.geno, pop1.inds, pop2.inds, SNPs, wHs)
    rm(Ht.dat, Ht.11, Ht.12, Ht.22, Ht.n, Ht)
    rm(Hs1.dat, Hs1.11, Hs1.12, Hs1.22, Hs1.n, Hs1)
    rm(Hs2.dat, Hs2.11, Hs2.12, Hs2.22, Hs2.n, Hs2)
    return(fst2)
  
}# End of function

# Apply function to data
cat("Calculate Fst\n")
fst <- llply(SNPs, function(x) if(length(x)>0) {
        func_Fst(data@gtdata,
                data@phdata[which(data@phdata$pop==1),]$id,
                data@phdata[which(data@phdata$pop==2),]$id,
                x)} else {c(0)}, .progress = "text",
    		.inform = TRUE, .parallel = FALSE, .paropts = NULL)

bed[,5] <- unlist(fst)
names(bed)[5] <- "Fst_mean"

# Write results to file
fname <- paste(args[8], "/Fst-chr", var_chr, "-win-", var_winL, ".txt", sep="")
write.table(bed, fname, sep="\t", row.names=T, col.names=T, append=F, quote=F)













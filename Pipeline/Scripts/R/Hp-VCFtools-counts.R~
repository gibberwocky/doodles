#!/usr/bin/env Rscript
args <- commandArgs(T)

## Generate allele frequency counts
## vcftools --vcf S2-07112014-OTH.vcf --out S2-07112014-OTH --counts
## vcftools --vcf S2-07112014-JFM.vcf --out S2-07112014-JFM --counts

## To call from terminal:
## ${RSCRIPTS}/Hp-VCFtools-counts.R ${arg1} ${arg2} ${arg3} ${arg4} ${arg5}

## arguments:
## [1] Directory containing input file 
## [2] Input file name (output from vcftools --count)
## [3] Window size in bp
## [4] Window step in bp
## [5] Tab-delimited file, 2 columns: chromosome name and length in bp (no header)
## Output file is ${arg1}/${arg2}-win-${arg3}.hp

## args <- c("/work/dwragg/Analysis/selection/hp", "S2-07112014-OTH.frq.count", "5000", "1000", "/save/dwragg/Apis/chrInfo")

# Set working directory
setwd(args[1])

# Read frequency counts data and remove SNPs with 0 allele data
dat <- read.table(args[2], sep="\t", row.names=NULL, header=T, strings=F)
names(dat) <- c("CHR", "POS", "N_ALLELES", "N_CHR", "A", "B")
dat <- dat[which(dat$N_CHR > 0),]

# Copy allele counts to new columns
dat[,7] <- as.numeric(gsub(".*:", "", dat$A))
dat[,8] <- as.numeric(gsub(".*:", "", dat$B))
names(dat)[7] <- "A_count"
names(dat)[8] <- "B_count"

# Read chromosome size data
chrInfo <- read.table(args[5], sep="\t", header=F, strings=F)


# ==============================================================================
# Generate windows for each chromosome
# ==============================================================================

bed <- NULL
for(i in chrInfo$V1){

	# Get length
	var_chr <- chrInfo[i,]$V1
	var_chrL <- chrInfo[i,]$V2

	# Specify window size
	var_winL <- as.numeric(args[3])
	var_overlap <- as.numeric(args[4])

	# Calculates intervals based on parameters
	start <- as.integer(seq(1,var_chrL, by=var_overlap))
	chr <- rep(var_chr, length(start))
	stop <- as.integer(start + var_winL)

	# Bind bed
	bed.tmp <- data.frame(chr, start, stop, stringsAsFactors=F)
	bed.tmp[nrow(bed.tmp),ncol(bed.tmp)] <- var_chrL
	bed <- rbind(bed, bed.tmp)

}


# ==============================================================================
# Calculate sums for major and minor alleles for SNPs in each window
# ==============================================================================

# Too slow to process all data at once, partition by chromosome
for(i in chrInfo$V1){

	# Subset data by chromosome
	tmp <- dat[which(dat$CHR==i),]
	tmp.bed <- bed[which(bed$chr==i),]

	# Identify SNPs per window
	nSNPs <- lapply(rownames(tmp.bed), function(x) 
		nrow(tmp[which(tmp$POS >= tmp.bed[x,]$start & 
			tmp$POS <= tmp.bed[x,]$stop),]) )

	# Sum Major
	sMaj <- lapply(rownames(tmp.bed), function(x) 
		sum(as.integer(do.call(pmax, tmp[which(tmp$POS >= tmp.bed[x,]$start & 
			tmp$POS <= tmp.bed[x,]$stop), 7:8],))) )

	# Sum Minor
	sMin <- lapply(rownames(tmp.bed), function(x) 
		sum(as.integer(do.call(pmin, tmp[which(tmp$POS >= tmp.bed[x,]$start & 
			tmp$POS <= tmp.bed[x,]$stop), 7:8],))) )

	# Update bed dataframe
	bed[rownames(tmp.bed),4] <- unlist(nSNPs)
	bed[rownames(tmp.bed),5] <- unlist(sMaj)
	bed[rownames(tmp.bed),6] <- unlist(sMin)
}

names(bed)[4] <- "nSNPs"
names(bed)[5] <- "sum_Maj"
names(bed)[6] <- "sum_Min"



# ==============================================================================
# Calculate pooled heterozygosity from major/minor allele frequencies
# ==============================================================================

# Calculate and record pooled heterozgosity = (2*aMj*aMn) / (aMj+aMn)^2
bed[,7] <- round(( 2 * bed[,5] * bed[,6] ) / ( bed[,5] + bed[,6] ) ^2, 3)
names(bed)[7] <- "Hp"

# Standardize Hp using mean (if normally distributed) or median (if skewed)
# Negative is sweep (positive selection)
# Positive is hypervariable (balancing selection)
var_Mean <- mean(bed$Hp, na.rm=T)
var_Median <- median(bed$Hp, na.rm=T)
var_Sdev <- sd(bed$Hp, na.rm=T)
bed[,8] <- round((bed[,7]-var_Mean)/var_Sdev, 3)
names(bed)[8] <- "Z_Hp"


# Write results to file
fname <- paste(args[1], "/", args[2], "-win-", var_winL, ".hp", sep="")
write.table(bed, fname, sep="\t", row.names=T, col.names=T, append=F, quote=F)




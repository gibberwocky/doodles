#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(IRanges)
library(stats)
library(pracma)
library(fdrtool)
library(mnormt)

# Read annotation data
genesA <- read.table("/home/dave/Copy/HoneyBee/Annotation/amel_OGSv3.2.gff3", sep="\t", header=F, strings=F, fill=T)
names(genesA) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
genesA <- genesA[which(genesA$type=="gene"),]
tmp.id <- gsub("ID=", "", genesA$attribute)    
table(genesB$gene%in%tmp.id)

genesB <- read.table("/home/dave/Copy/HoneyBee/Annotation/Genes-BioMartExport_incStrand.txt", sep="\t", header=T, strings=F, fill=T)
names(genesB) <- c("gene", "chr", "start", "end")

ir.genes <- GRanges(seqnames=genes$chr, ranges=IRanges(start=genes$start, end=genes$end))

# Set path
setwd("/home/dave/Documents/SeqApiPop")
#windows <- c("1kb", "2kb", "3kb", "4kb", "5kb", "10kb", "20kb", "50kb")
windows <- c("3kb", "5kb", "10kb", "20kb")

# Read all significant results in
results <- list()
Fst <- list()
Hp <- list()
Hl <- list()
for(i in 1:length(windows)){
  fname <- paste("Results/", windows[i], "/significant.bed", sep="")
  results[[i]] <- read.table(fname, sep="\t", header=F, strings=F)
  names(results[[i]]) <- c("chr", "start", "end", "Val", "method")
  Fst[[i]] <- results[[i]][which(results[[i]]$method=="Fst"),]
  Hp[[i]] <- results[[i]][grep("Hp", results[[i]]$method),]
  Hl[[i]] <- results[[i]][grep("Local", results[[i]]$method),]
}

# Filter results
pop <- ""
tag <- "Fst"
qval <- 0
tmp.ir <- NULL
tmp.irs <- NULL
tmp.genes <- NULL
tmp.results <- NULL
for(i in 1:length(windows)){
  tmp <- results[[i]][grep(pop, results[[i]]$method),]
  tmp <- tmp[grep(tag, tmp$method),]
  tmp <- tmp[which(tmp$Val >= quantile(tmp$Val, qval)),]
  
  # Convert to genomic range
  tmp.ir <- reduce(GRanges(seqnames=tmp$chr, ranges=IRanges(start=tmp$start, end=tmp$end)))

  # Within the given window-size dataset, record number of windows merged and mean value
  tmp.dat <- data.frame(val=numeric(1), num=numeric(1), stringsAsFactors=F)
  for(n in 1:length(tmp.ir)){
    a <- tmp[which(tmp$chr==as.character(seqnames(tmp.ir[n])) &
              tmp$start >= tmp.ir@ranges@start[n] &
              tmp$end <= (tmp.ir@ranges@start[n] + tmp.ir@ranges@width[n]-1)),]$Val
    tmp.dat[n,1] <- round(mean(abs(a), na.rm=T), 3)
    tmp.dat[n,2] <- length(a)
  }
  mcols(tmp.ir) <- tmp.dat
  
  # Filter out intervals supported by only 1 window
  tmp.ir <- tmp.ir[which(tmp.ir$num>1),]

  # Identify overlapping genes
  tmp.genes <- rbind(tmp.genes, genes[findOverlaps(tmp.ir, ir.genes)@subjectHits,])
  
  # Record reduced window details
  tmp <- data.frame(chr=seqnames(tmp.ir), start=tmp.ir@ranges@start, end=tmp.ir@ranges@start + tmp.ir@ranges@width - 1, 
                    val=tmp.ir$val, num=tmp.ir$num, stringsAsFactors=F)
  tmp.irs <- rbind(tmp.irs, tmp)
}

# Identify genes that occurred in more than one window size (in other words, are duplicated)
tmp <- tmp.genes[duplicated(tmp.genes$gene),]
# Get unique gene IDs from the duplicates
tmp <- tmp[!duplicated(tmp$gene),]

# Write genes to file
fname <- paste(tag, ".genes", sep="")
write.table(tmp, fname, append=F, quote=F, row.names=F, col.names=T, sep="\t")

# Reduce intervals across all window size datasets
tmp.rirs <- reduce(GRanges(seqnames=tmp.irs$chr, ranges=IRanges(start=tmp.irs$start, end=tmp.irs$end)))
tmp.dat <- data.frame(val=numeric(1), num=numeric(1), stringsAsFactors=F)
for(n in 1:length(tmp.rirs)){
  a <- tmp.irs[which(tmp.irs$chr==as.character(seqnames(tmp.rirs[n])) &
                   tmp.irs$start >= tmp.rirs@ranges@start[n] &
                   tmp.irs$end <= (tmp.rirs@ranges@start[n] + tmp.rirs@ranges@width[n]-1)),]$val
  tmp.dat[n,1] <- round(mean(abs(a), na.rm=T), 3)
  tmp.dat[n,2] <- length(a)
}
mcols(tmp.rirs) <- tmp.dat

# Write intervals to file
fname <- paste(tag, ".irs", sep="")
tmp.dat <- data.frame(chr=seqnames(tmp.rirs), start=tmp.rirs@ranges@start, end=(tmp.rirs@ranges@start+tmp.rirs@ranges@width-1),
                      val=tmp.rirs$val, num=tmp.rirs$num, stringsAsFactors=F)
write.table(tmp.dat, fname, append=F, quote=F, row.names=F, col.names=T, sep="\t")





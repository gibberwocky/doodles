# Load haplotype plotting function
source("/home/dave/Copy/HoneyBee/Notes/Pipeline/R/Function-HaploPlot.R")

# Read genome annotation data
genes <- read.table("/home/dave/Copy/HoneyBee/Annotation/Apis_mellifera.GCA_000002195.1.24.gtf", sep="\t", header=F, strings=F, fill=T)
names(genes) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Read sample lists
cds.path <- "/home/dave/Copy/HoneyBee/Analyses/selection/CDS"
JFM <- unlist(read.table(paste(cds.path, "/HAPS/JFM_samples_fixed.list", sep=""), header=F, stringsAsFactors=F))
OTH <- unlist(read.table(paste(cds.path, "/HAPS/OTH_samples_fixed.list", sep=""), header=F, stringsAsFactors=F))

# JFM non-synonymous hap data
JFM.N.pos <- read.table(paste(cds.path, "/HAPS/JFM-N.snps.012.pos", sep=""), header=F, stringsAsFactors=F)
names(JFM.N.pos) <- c("chr", "pos")
JFM.N.hap <- data.frame(t(read.table(paste(cds.path, "/HAPS/JFM-N.snps.012", sep=""), header=F, stringsAsFactors=F))[-1,])
names(JFM.N.hap) <- unlist(read.delim(paste(cds.path, "/HAPS/JFM-N.snps.012.indv", sep=""), header=F, sep="\n", stringsAsFactors=F))
JFM.N.hap <- cbind(JFM.N.pos, JFM.N.hap)
rm(JFM.N.pos)

# JFM synonymous hap data
JFM.S.pos <- read.table(paste(cds.path, "/HAPS/JFM-S.snps.012.pos", sep=""), header=F, stringsAsFactors=F)
names(JFM.S.pos) <- c("chr", "pos")
JFM.S.hap <- data.frame(t(read.table(paste(cds.path, "/HAPS/JFM-S.snps.012", sep=""), header=F, stringsAsFactors=F))[-1,])
names(JFM.S.hap) <- unlist(read.delim(paste(cds.path, "/HAPS/JFM-S.snps.012.indv", sep=""), header=F, sep="\n", stringsAsFactors=F))
JFM.S.hap <- cbind(JFM.S.pos, JFM.S.hap)
rm(JFM.S.pos)

# OTH non-synonymous hap data
OTH.N.pos <- read.table(paste(cds.path, "/HAPS/OTH-N.snps.012.pos", sep=""), header=F, stringsAsFactors=F)
names(OTH.N.pos) <- c("chr", "pos")
OTH.N.hap <- data.frame(t(read.table(paste(cds.path, "/HAPS/OTH-N.snps.012", sep=""), header=F, stringsAsFactors=F))[-1,])
names(OTH.N.hap) <- unlist(read.delim(paste(cds.path, "/HAPS/OTH-N.snps.012.indv", sep=""), header=F, sep="\n", stringsAsFactors=F))
OTH.N.hap <- cbind(OTH.N.pos, OTH.N.hap)
rm(OTH.N.pos)

# OTH synonymous hap data
OTH.S.pos <- read.table(paste(cds.path, "/HAPS/OTH-S.snps.012.pos", sep=""), header=F, stringsAsFactors=F)
names(OTH.S.pos) <- c("chr", "pos")
OTH.S.hap <- data.frame(t(read.table(paste(cds.path, "/HAPS/OTH-S.snps.012", sep=""), header=F, stringsAsFactors=F))[-1,])
names(OTH.S.hap) <- unlist(read.delim(paste(cds.path, "/HAPS/OTH-S.snps.012.indv", sep=""), header=F, sep="\n", stringsAsFactors=F))
OTH.S.hap <- cbind(OTH.S.pos, OTH.S.hap)
rm(OTH.S.pos)

# Read genesets into a table
dat.path <- "/home/dave/Documents/SeqApiPop/GeneSets"
for(i in 1:6){
  if(i==1) fname <- paste(dat.path, "/DivergentGenes", sep="")
  if(i==2) fname <- paste(dat.path, "/HoneyGenes", sep="")
  if(i==3) fname <- paste(dat.path, "/JellyGenes", sep="")
  if(i==4) fname <- paste(dat.path, "/DivergentGenesNoGO", sep="")
  if(i==5) fname <- paste(dat.path, "/HoneyGenesNoGO", sep="")
  if(i==6) fname <- paste(dat.path, "/JellyGenesNoGO", sep="")
  tmp <- read.table(fname, header=F, sep="\t", stringsAsFactors=F)
  if(i==1) hits <- tmp
  if(i>1) hits <- rbind(hits, tmp)
}
hits <- hits[!duplicated(hits[,1]),]
for(i in 1:nrow(hits)){
  gene <- genes[grep(hits[i,1], genes$attribute),]
  hits[i,2] <- gene$chr[1]
  hits[i,3] <- min(gene$start)
  hits[i,4] <- max(gene$end)
}
names(hits) <- c("gene", "chr", "bp1", "bp2")



# output path
plot.path <- "/home/dave/Documents/SeqApiPop/Plots/Genes"

# Loops through each chromosome and produces plot for each hit on the chromosome from the CDS analysis
for(chromosome in unique(hits$chr)){
  
  # Read map file
  fname <- paste("/home/dave/Copy/HoneyBee/Analyses/haploblocks/current/CHR/", chromosome, ".map", sep="")
  map.haps <- read.table(fname, header=F, stringsAsFactors=F)
  names(map.haps)[4] <- "start"
  
  # Read JFM haplotype files
  fname <- paste("/home/dave/Copy/HoneyBee/Analyses/haploblocks/current/JFM/JFM-", chromosome, ".hap", sep="")
  jfm.haps <- read.table(fname, header=F, stringsAsFactors=F)
  names(jfm.haps) <- map.haps$start
  
  # Read OTH haplotype files
  fname <- paste("/home/dave/Copy/HoneyBee/Analyses/haploblocks/current/OTH/OTH-", chromosome, ".hap", sep="")
  oth.haps <- read.table(fname, header=F, stringsAsFactors=F)
  names(oth.haps) <- map.haps$start
  
  # Loop through each hit on chromosome
  dat.hits <- hits[which(hits$chr==chromosome),]
  for(hit in 1:nrow(dat.hits)){
    
    # Get annotation details for gene
    target <- dat.hits[hit,1]
    gene <- genes[grep(target, genes$attribute),]
    chr <- gene$chr[1]
    bp1 <- min(gene$start)
    bp2 <- max(gene$end)
    
    # Subset non-synonymous and synonmous details to target region
    hap1.N <- JFM.N.hap[which(JFM.N.hap$chr == chr & JFM.N.hap$pos >= bp1 & JFM.N.hap$pos <= bp2),]
    hap1.S <- JFM.S.hap[which(JFM.S.hap$chr == chr & JFM.S.hap$pos >= bp1 & JFM.S.hap$pos <= bp2),]
    hap2.N <- OTH.N.hap[which(OTH.N.hap$chr == chr & OTH.N.hap$pos >= bp1 & OTH.N.hap$pos <= bp2),]
    hap2.S <- OTH.S.hap[which(OTH.S.hap$chr == chr & OTH.S.hap$pos >= bp1 & OTH.S.hap$pos <= bp2),]  
    # Merge the data from both populations into single N and S data frames
    hap.N <- rbind(hap1.N, hap2.N)
    hap.N <- hap.N[!duplicated(hap.N),]
    hap.S <- rbind(hap1.S, hap2.S)
    hap.S <- hap.S[!duplicated(hap.S),]    
    
    # Plot gene and alleles
    fname <- paste(plot.path, "/", chromosome, "-", bp1, "-", bp2, "_", target, ".png", sep="")
    p <- haploplot(dat.hits[hit,1], gene, jfm.haps, oth.haps, hap.N, hap.S, bp1, bp2, 500, "HN", "RJ", JFM, OTH,
                  fname, "png", units="mm", width=169, height=169, res=300)  
    if(p==0) cat(chromosome, ":", dat.hits[hit,]$bp1, "-", dat.hits[hit,]$bp2, " returned no plot", "\n")
    # Plots can be empty if the extension from the gene is not large, and therefore not proximal to Fst hit
    
  }
  
}

  
  
  
  
  
  
  
  


#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(IRanges)
library(plotrix)

# Set path
outpath <- "/home/dave/Documents/SeqApiPop"
dir.create(file.path(outpath, paste("Results/Coverage", sep="")), recursive=T)
win_size <- "5kb"
snp_min <- 1


setwd("/home/dave/Copy/HoneyBee/Analyses/selection")

# Depth of coverage JFM
fname.dp.jfm <- paste("DP/", win_size, "/JFM.dp", sep="")
dat.dp.jfm <- read.table(fname.dp.jfm, header=T, stringsAsFactors=F)
names(dat.dp.jfm)[2] <- "start"
names(dat.dp.jfm)[5] <- "Val"
dat.dp.jfm <- dat.dp.jfm[which(dat.dp.jfm$Val < (mean(dat.dp.jfm$Val) + 3*sd(dat.dp.jfm$Val))),]

# Depth of coverage OTH
fname.dp.oth <- paste("DP/", win_size, "/OTH.dp", sep="")
dat.dp.oth <- read.table(fname.dp.oth, header=T, stringsAsFactors=F)
names(dat.dp.oth)[2] <- "start"
names(dat.dp.oth)[5] <- "Val"
dat.dp.oth <- dat.dp.oth[which(dat.dp.oth$Val < (mean(dat.dp.oth$Val) + 3*sd(dat.dp.oth$Val))),]


# VCF-Nucleotide Diversity JFM
fname.pi.jfm <- paste("Nuc-Div/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle-JFM.windowed.pi", sep="")
dat.pi.jfm <- read.table(fname.pi.jfm, header=T, stringsAsFactors=F)
names(dat.pi.jfm)[1] <- "chr"
names(dat.pi.jfm)[2] <- "start"
names(dat.pi.jfm)[5] <- "Val"

# VCF-Nucleotide Diversity OTH
fname.pi.oth <- paste("Nuc-Div/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle-OTH.windowed.pi", sep="")
dat.pi.oth <- read.table(fname.pi.oth, header=T, stringsAsFactors=F)
names(dat.pi.oth)[1] <- "chr"
names(dat.pi.oth)[2] <- "start"
names(dat.pi.oth)[5] <- "Val"


# ================================================================================================
# Plot of nucleotide diversity versus depth of coverage
# ================================================================================================
for(chromosome in seq(1,16)) {

  dat.dp.jfm.tmp <- dat.dp.jfm[which(dat.dp.jfm$chr == chromosome),]
  dat.dp.oth.tmp <- dat.dp.oth[which(dat.dp.oth$chr == chromosome),]
  dat.pi.jfm.tmp <- dat.pi.jfm[which(dat.pi.jfm$chr == chromosome),]
  dat.pi.oth.tmp <- dat.pi.oth[which(dat.pi.oth$chr == chromosome),]
  
  tmp.jfm <- merge(dat.dp.jfm.tmp, dat.pi.jfm.tmp, by="start")
  tmp.oth <- merge(dat.dp.oth.tmp, dat.pi.oth.tmp, by="start")
  
  # Plot parameters
  setwd(outpath)
  pdfname <- paste("Results/Coverage/", win_size, "/", chromosome, "-Pi-DP.png", sep="")
  png(pdfname, height=9, width=6, units="in", res=300)
  
  # Set-up
  bkgr="#FFFFFF"
  frgr="#000000"
  hicov=rgb(1.0,0.0,0.0)
  locov=rgb(1.0,0.0,0.0)
  par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
      cex=0.3, cex.sub=0.7, cex.main=0.7, cex.axis=0.7, cex.lab=0.7, lwd=1)
  vrow <- 2
  vcol <- 1
  par(mfrow=c(vrow,vcol), oma = c(2,2,1,1), mar = c(6,4,2,2))
  pt.size <- 0.2
  
  # Plots
  for(i in seq(1,2)){
  if(i==1) { tmp <- tmp.jfm; main_txt = "Honey production" }
  if(i==2) { tmp <- tmp.oth; main_txt = "Royal jelly production" }
  plot(NULL, xlab = "DP", ylab = expression(pi), main = main_txt, 
       xlim = c(min(tmp.jfm$Val.x, tmp.oth$Val.x), max(tmp.jfm$Val.x, tmp.oth$Val.x)), 
       ylim = c(min(tmp.jfm$Val.y, tmp.oth$Val.y), max(tmp.jfm$Val.y, tmp.oth$Val.y)))
  with(tmp, points(x=tmp$Val.x, y=tmp$Val.y, col = rgb(0,0,0.6), pch=19, cex=pt.size))
  draw.ellipse(x = mean(tmp$Val.x), y = mean(tmp$Val.y), a = 2*sd(tmp$Val.x), b = 2*sd(tmp$Val.y),
               lwd=3, border=rgb(1,0.3,0.3), lty=1)
  }  
  title(xlab = paste("Chromosome ", chromosome, sep=""), line=-0.5, outer=T)
  
  dev.off()
}



# ================================================================================================
# Plot of depth of coverage
# ================================================================================================
for(chromosome in seq(1,16)) {
  
  dat.dp.jfm.tmp <- dat.dp.jfm[which(dat.dp.jfm$chr == chromosome),]
  dat.dp.oth.tmp <- dat.dp.oth[which(dat.dp.oth$chr == chromosome),]  
  
  # Plot parameters
  setwd(outpath)
  pdfname <- paste("Results/Coverage/", win_size, "/", chromosome, "-DP.png", sep="")
  png(pdfname, height=6, width=6, units="in", res=300)
  
  # Set-up
  bkgr="#FFFFFF"
  frgr="#000000"
  hicov=rgb(1.0,0.0,0.0)
  locov=rgb(1.0,0.0,0.0)
  par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
      cex=0.3, cex.sub=0.7, cex.main=0.7, cex.axis=0.7, cex.lab=0.7, lwd=1)
  vrow <- 2
  vcol <- 1
  par(mfrow=c(vrow,vcol), oma = c(5.5,1.5,1,1), mar = c(1.5,3,1,0))
  xmin <- 0
  xmax <- max(dat.dp.jfm.tmp$end, dat.dp.oth.tmp$end)
  ymin <- 0
  ymax <- max(dat.dp.jfm$Val, dat.dp.oth$Val)+100
  pt.size <- 0.2  
  
  for(i in seq(1,2)){
    if(i==1) tmp <- dat.dp.jfm.tmp
    if(i==2) tmp <- dat.dp.oth.tmp
    
    # Plot
    plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
         xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = NA, 
         ylab = "", las = 1, pch = 19)
    mtext("DP", 2, line=3, cex=0.7)
    
    nchr <- 17
    col <- rep(rgb(0.85,0.85,0.85), nchr)
    col[seq(1,nchr,2)] <- rgb(0.8,0.8,0.8)
    icol <- 1    
  
    if(i==1) { title(ylab = "", outer=F, line= 2); title (main="Depth of coverage (HN)", line = 0.4) }
    if(i==2) { title(ylab = "", outer=F, line= 2); title (main="Depth of coverage (RJ)", line = 0.4) }
    
    # Significance thresholds and mean indicator
    with(tmp, points(x=start, y=Val, col = col[icol], pch=19, cex=pt.size))
    with(tmp[which(tmp$Val < 10),], points(x=start, y=Val, col = locov, pch=19, cex=pt.size))
    with(tmp[which(tmp$Val > (3*mean(tmp$Val))),], points(x=start, y=Val, col = hicov, pch=19, cex=pt.size))
    abline(h=round(mean(tmp$Val, na.rm=T),3), col=rgb(0.3,0.3,0.3), lty=2, lwd=1)    
    abline(h=10, col=locov, lty=2, lwd=1)    
    abline(h=3*mean(tmp$Val), col=hicov, lty=2, lwd=1)    
    
    box()
   
  }
  
  title(xlab = paste("Chromosome ", chromosome, " position", sep=""), outer=T)
  axis(1, col=frgr, at=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))],
       labels=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))]) 
  
  dev.off()

}






# ================================================================================================
# Summary statistics
# ================================================================================================

# Correlation tests
# DP~PI; JFM r = 0.0349  OTH r = 0.0347
# GC~PI; JFM r = 0.0233  OTH r = 0.0223
# GC~DP; JFM r = 0.0833  OTH r = 0.0732

# Mean (sd) values
# DP; JFM 217.48 (65.15)  OTH 189.8 (69.48)
# GC; JFM 0.309 (0.12)  OTH 0.309 (0.12)
# Pi; JFM 0.0038 (0.0022)  OTH 0.0035 (0.0022)

# Hi (.99 quantile) and Lo (.1 quantile) diversity regions in each population
JFM.lopi <- dat.pi.jfm[which(dat.pi.jfm$Val < quantile(dat.pi.jfm$Val, .01)),]
JFM.lopi <- reduce(GRanges(seqnames=JFM.lopi$chr, ranges=IRanges(start=JFM.lopi$start, end=JFM.lopi$BIN_END))) # 613 intervals
write.table(paste(seqnames(JFM.lopi), start(JFM.lopi), end(JFM.lopi), sep=":"), file=paste("Pi/", win_size, "/JFM.lopi", sep=""), 
            append=F, quote=F, row.names=F, col.names=F)
JFM.hipi <- dat.pi.jfm[which(dat.pi.jfm$Val > quantile(dat.pi.jfm$Val, .99)),]
JFM.hipi <- reduce(GRanges(seqnames=JFM.hipi$chr, ranges=IRanges(start=JFM.hipi$start, end=JFM.hipi$BIN_END))) # 652 intervals
write.table(paste(seqnames(JFM.hipi), start(JFM.hipi), end(JFM.hipi), sep=":"), file=paste("Pi/", win_size, "/JFM.hipi", sep=""), 
            append=F, quote=F, row.names=F, col.names=F)
OTH.lopi <- dat.pi.oth[which(dat.pi.oth$Val < quantile(dat.pi.oth$Val, .01)),]
OTH.lopi <- reduce(GRanges(seqnames=OTH.lopi$chr, ranges=IRanges(start=OTH.lopi$start, end=OTH.lopi$BIN_END))) # 742 intervals
write.table(paste(seqnames(OTH.lopi), start(OTH.lopi), end(OTH.lopi), sep=":"), file=paste("Pi/", win_size, "/OTH.lopi", sep=""), 
            append=F, quote=F, row.names=F, col.names=F)
OTH.hipi <- dat.pi.oth[which(dat.pi.oth$Val > quantile(dat.pi.oth$Val, .99)),]
OTH.hipi <- reduce(GRanges(seqnames=OTH.hipi$chr, ranges=IRanges(start=OTH.hipi$start, end=OTH.hipi$BIN_END))) # 586 intervals
write.table(paste(seqnames(OTH.hipi), start(OTH.hipi), end(OTH.hipi), sep=":"), file=paste("Pi/", win_size, "/OTH.hipi", sep=""), 
            append=F, quote=F, row.names=F, col.names=F)





# ================================================================================================
# Plot of gene count versus nucleotide diversity for all window genome-wide
# ================================================================================================

# Read annotation data
genes <- read.table("/home/dave/Copy/HoneyBee/Annotation/Genes-BioMartExport.txt", sep="\t", header=T, strings=F, fill=T)
names(genes) <- c("gene", "chr", "start", "end")

# Plot parameters
setwd(outpath)
pdfname <- paste("Results/Coverage/", win_size, "/GeneCount-Pi.png", sep="")
png(pdfname, height=6, width=6, units="in", res=300)

# Set-up
bkgr="#FFFFFF"
frgr="#000000"
hicov=rgb(1.0,0.0,0.0)
locov=rgb(1.0,0.0,0.0)
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
    cex=0.3, cex.sub=0.7, cex.main=0.7, cex.axis=0.7, cex.lab=0.7, lwd=1)
vrow <- 2
vcol <- 1
par(mfrow=c(vrow,vcol), oma = c(2,2,2,1), mar = c(3,3,3,2))

main_txt <- ""
pt.size <- 0.2

# Plots
for(i in seq(1,2)){
  if(i==1) { dat <- dat.pi.jfm; main_txt = "Window gene density versus nucleotide diversity (HN)" }
  if(i==2) { dat <- dat.pi.oth; main_txt = "Window gene density versus nucleotide diversity (RJ)" }
  ir1 <- GRanges(seqnames=genes$chr, ranges=IRanges(start=genes$start, end=genes$end))
  ir2 <- GRanges(seqnames=dat$chr, ranges=IRanges(start=dat$start, end=dat$BIN_END), pi=dat$Val)
  dat[,6] <- countOverlaps(ir2, ir1)
  names(dat)[6] <- "Genes"
  plot(NULL, ylab = "Gene count", xlab = expression(pi), main = main_txt, 
     xlim = c(min(dat$Val), max(dat$Val)), ylim = c(min(dat$Genes), max(dat$Genes)))
  with(dat, points(x=dat$Val, y=dat$Genes, col = rgb(0,0,0.6), pch=19, cex=pt.size))
  title(xlab = expression(pi), line=2, outer=F)
  title(ylab = "Gene count", line=2, outer=F)
}

dev.off()



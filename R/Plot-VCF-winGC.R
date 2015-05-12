library(plotrix)

# Set path
outpath <- "/home/dave/Documents/SeqApiPop"
dir.create(file.path(outpath, paste("Results/GCcontent", sep="")), recursive=T)
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

# GC content JFM
fname.gc.jfm <- paste("GC/", win_size, "/JFM.gc", sep="")
dat.gc.jfm <- read.table(fname.gc.jfm, header=F, stringsAsFactors=F)
names(dat.gc.jfm) <- c("chr", "start", "end", "GC")

# GC content OTH
fname.gc.oth <- paste("GC/", win_size, "/OTH.gc", sep="")
dat.gc.oth <- read.table(fname.gc.oth, header=F, stringsAsFactors=F)
names(dat.gc.oth) <- c("chr", "start", "end", "GC")

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
# Plot of GC content versus nucleotide diversity
# ================================================================================================
for(chromosome in seq(1,16)) {

  dat.gc.jfm.tmp <- dat.gc.jfm[which(dat.gc.jfm$chr == chromosome),]
  dat.gc.oth.tmp <- dat.gc.oth[which(dat.gc.oth$chr == chromosome),]
  dat.pi.jfm.tmp <- dat.pi.jfm[which(dat.pi.jfm$chr == chromosome),]  
  dat.pi.oth.tmp <- dat.pi.oth[which(dat.pi.oth$chr == chromosome),]
  
  tmp.jfm <- merge(dat.gc.jfm.tmp, dat.pi.jfm.tmp, by="start")
  tmp.oth <- merge(dat.gc.oth.tmp, dat.pi.oth.tmp, by="start")
  
  # Plot parameters
  setwd(outpath)
  pdfname <- paste("Results/GCcontent/", win_size, "/", chromosome, "-Pi-GC.png", sep="")
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
    plot(NULL, xlab = "GC", ylab = expression(pi), main = main_txt, 
         xlim = c(min(tmp.jfm$GC, tmp.oth$GC), max(tmp.jfm$GC, tmp.oth$GC)), 
         ylim = c(min(tmp.jfm$Val, tmp.oth$Val), max(tmp.jfm$Val, tmp.oth$Val)))
    with(tmp, points(x=tmp$GC, y=tmp$Val, col = rgb(0,0,0.6), pch=19, cex=pt.size))
    draw.ellipse(x = mean(tmp$GC), y = mean(tmp$Val), a = 2*sd(tmp$GC), b = 2*sd(tmp$Val),
                 lwd=3, border=rgb(1,0.3,0.3), lty=1)
  }  
  title(xlab = paste("Chromosome ", chromosome, sep=""), line=-0.5, outer=T)
  
  dev.off()
}





# ================================================================================================
# Plot of GC content versus depth of coverage
# ================================================================================================
for(chromosome in seq(1,16)) {
  
  dat.dp.jfm.tmp <- dat.dp.jfm[which(dat.dp.jfm$chr == chromosome),]
  dat.dp.oth.tmp <- dat.dp.oth[which(dat.dp.oth$chr == chromosome),]
  dat.gc.jfm.tmp <- dat.gc.jfm[which(dat.gc.jfm$chr == chromosome),]
  dat.gc.oth.tmp <- dat.gc.oth[which(dat.gc.oth$chr == chromosome),]

  setwd("/home/dave/Copy/HoneyBee/Analyses/selection")
    
  tmp.jfm <- merge(dat.gc.jfm.tmp, dat.dp.jfm.tmp, by="start")
  tmp.oth <- merge(dat.gc.oth.tmp, dat.dp.oth.tmp, by="start")
  
  # Plot parameters
  setwd(outpath)
  pdfname <- paste("Results/GCcontent/", win_size, "/", chromosome, "-GC-DP.png", sep="")
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
    plot(NULL, xlab = "GC", ylab = "DP", main = main_txt, 
         xlim = c(min(tmp.jfm$GC, tmp.oth$GC), max(tmp.jfm$GC, tmp.oth$GC)), 
         ylim = c(min(tmp.jfm$Val, tmp.oth$Val), max(tmp.jfm$Val, tmp.oth$Val)))
    with(tmp, points(x=tmp$GC, y=tmp$Val, col = rgb(0,0,0.6), pch=19, cex=pt.size))
    draw.ellipse(x = mean(tmp$GC), y = mean(tmp$Val), a = 2*sd(tmp$GC), b = 2*sd(tmp$Val),
                 lwd=3, border=rgb(1,0.3,0.3), lty=1)
  }  
  title(xlab = paste("Chromosome ", chromosome, sep=""), line=-0.5, outer=T)
  
  dev.off()
}




# ================================================================================================
# Plot of GC content
# ================================================================================================
for(chromosome in seq(1,16)) {
  
  dat.gc.jfm.tmp <- dat.gc.jfm[which(dat.gc.jfm$chr == chromosome),]
  dat.gc.oth.tmp <- dat.gc.oth[which(dat.gc.oth$chr == chromosome),]  
  
  
  # Plot parameters
  setwd(outpath)
  pdfname <- paste("Results/GCcontent/", win_size, "/", chromosome, "-GC.png", sep="")
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
  par(mfrow=c(vrow,vcol), oma = c(5.5,1.5,1,1), mar = c(1.5,3,1,1))
  xmin <- 0
  xmax <- max(dat.gc.jfm.tmp$end, dat.gc.oth.tmp$end)
  ymin <- 0
  ymax <- max(dat.gc.jfm$GC, dat.gc.oth$GC) + 0.1
  pt.size <- 0.2
  
  for(i in seq(1,2)){
    if(i==1) { tmp <- dat.gc.jfm.tmp; tmp2 <- dat.gc.oth.tmp }
    if(i==2) { tmp <- dat.gc.oth.tmp; tmp2 <- dat.gc.jfm.tmp }
    tmp.d <- tmp[which(tmp$GC != tmp2$GC),]
    
    # Plot
    par(new=F)
    ymin <- 0
    ymax <- max(dat.gc.jfm$GC, dat.gc.oth$GC) + 0.1
    plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
         xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = NA, 
         ylab = "", las = 1, pch = 19)
    mtext("GC", 2, line=3, cex=0.7)
    
    nchr <-17
    col <- rep(rgb(0.85,0.85,0.85), nchr)
    col[seq(1,nchr,2)] <- rgb(0.8,0.8,0.8)
    icol <- 1    
  
    if(i==1) { title(ylab = "", outer=F, line= 2); title (main=paste("GC content in ", win_size, " windows (HN)", sep=""), line = 0.4) }
    if(i==2) { title(ylab = "", outer=F, line= 2); title (main=paste("GC content in ", win_size, " windows (RJ)", sep=""), line = 0.4) }
    
    # Significance thresholds and mean indicator
    with(tmp, points(x=start, y=GC, col = col[icol], pch=19, cex=pt.size))
    with(tmp[which(tmp$GC < (mean(tmp$GC) - 2*(sd(tmp$GC)))),], points(x=start, y=GC, col = locov, pch=19, cex=pt.size))
    with(tmp[which(tmp$GC > (mean(tmp$GC) + 2*(sd(tmp$GC)))),], points(x=start, y=GC, col = hicov, pch=19, cex=pt.size))
    with(tmp[which(tmp$GC != tmp2$GC),], points(x=start, y=GC, col=rgb(0.65,0.65,0.65), pch=19, cex=pt.size))
    with(tmp[which( (tmp$GC - tmp2$GC) > 0.005),], points(x=start, y=GC, col=rgb(0,0,0.8), pch=19, cex=pt.size))
    abline(h=round(mean(tmp$GC, na.rm=T),3), col=rgb(0.3,0.3,0.3), lty=2, lwd=1)    
    abline(h=mean(tmp$GC) - 2*(sd(tmp$GC)), col=locov, lty=2, lwd=1)    
    abline(h=mean(tmp$GC) + 2*(sd(tmp$GC)), col=hicov, lty=2, lwd=1)    
    
    box()
   
  }
  
  title(xlab = paste("Chromosome ", chromosome, " position", sep=""), outer=T)
  axis(1, col=frgr, at=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))],
       labels=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))]) 
  
  dev.off()

}








# ================================================================================================
# Plot of gene count versus nucleotide diversity for all window genome-wide
# ================================================================================================

# Read annotation data
genes <- read.table("/home/dave/Copy/HoneyBee/Annotation/Genes-BioMartExport.txt", sep="\t", header=T, strings=F, fill=T)
names(genes) <- c("gene", "chr", "start", "end")

# Plot parameters
setwd(outpath)
pdfname <- paste("Results/GCcontent/", win_size, "/GeneCount-GC.png", sep="")
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
pt.size <- 0.2

# Plots
for(i in seq(1,2)){
  if(i==1) { dat <- dat.gc.jfm; main_txt = "Window gene density versus GC content (HN)" }
  if(i==2) { dat <- dat.gc.oth; main_txt = "Window gene density versus GC content (RJ)" }
  ir1 <- GRanges(seqnames=genes$chr, ranges=IRanges(start=genes$start, end=genes$end))
  ir2 <- GRanges(seqnames=dat$chr, ranges=IRanges(start=dat$start, end=dat$end), gc=dat$GC)
  dat[,5] <- countOverlaps(ir2, ir1)
  names(dat)[5] <- "Genes"
  plot(NULL, ylab = "Gene count", xlab = "GC content", main = main_txt, 
       xlim = c(min(dat$GC), max(dat$GC)), ylim = c(min(dat$Genes), max(dat$Genes)))
  with(dat, points(x=dat$GC, y=dat$Genes, col = rgb(0,0,0.6), pch=19, cex=pt.size))
  title(xlab = "GC content", line=2, outer=F)
  title(ylab = "Gene count", line=2, outer=F)
}

dev.off()

















#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(IRanges)
library(stats)
library(pracma)
library(fdrtool)

#fsttmp <- ( dat.fst$Val - mean(dat.fst$Val, na.rm=T) ) / sd(dat.fst$Val, na.rm=T)
#test <- fdrtool(fsttmp)

# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/selection")

outpath <- "/home/dave/Documents/SeqApiPop"


for(n in 1:3){
  
  if(n==1) pdfname <- paste(outpath, "/Results/hp-histograms.png", sep="")
  if(n==2) pdfname <- paste(outpath, "/Results/hp-histograms-JFM.png", sep="")
  if(n==3) pdfname <- paste(outpath, "/Results/hp-histograms-OTH.png", sep="")
  
  png(pdfname, height=6, width=10, units="in", res=300)
  ZUpper <- 4 # Z=3.719 =~ Q=0.0001 =~ 1 in 9999
  ZLower <- 3.719 # Z=3 =~ Q=0.001349 =~ 1 in 741)
  pcolz1=rgb(0.0,0.0,1.0)
  pcolz2=rgb(0.8,0.8,1.0)
  vrow <- 2
  vcol <- 3
  par(mfrow=c(vrow,vcol), oma = c(1,1,1,1), mar = c(2,2,2,2), cex=0.7, cex.axis=0.7)
  
  for(i in 1:6){
    if(i==1) { win_size <- "1kb"; size <- 1000; maintxt <- "A" }
    if(i==2) { win_size <- "5kb"; size <- 5000; maintxt <- "B" }
    if(i==3) { win_size <- "50kb"; size <- 50000; maintxt <- "C" }
    if(i==4) { win_size <- "1kb"; size <- 1000; maintxt <- "A" }
    if(i==5) { win_size <- "5kb"; size <- 5000; maintxt <- "B" }
    if(i==6) { win_size <- "50kb"; size <- 50000; maintxt <- "C" }
    
    # Pooled heterozygosity combined data
    if(n==1) fname.hp <- paste("Hp/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle.vcf.hp", sep="")
    if(n==2) fname.hp <- paste("Hp/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle-JFM.vcf.hp", sep="")
    if(n==3) fname.hp <- paste("Hp/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle-OTH.vcf.hp", sep="")
    dat.hp <- read.table(fname.hp, header=T, stringsAsFactors=F)
    dat.hp <- dat.hp[which(dat.hp$chr<17),]
    dat.hp <- dat.hp[which(dat.hp$n_snps>10),]
    dat.hp[,11] <- ( dat.hp$Hp - mean(dat.hp$Hp, na.rm=T) ) / sd(dat.hp$Hp, na.rm=T)
    names(dat.hp)[11] <- "ZHp"
    dat <- dat.hp
      
    if(i<4){
      xlim = c(quantile(dat$Hp,0), quantile(dat$Hp,1))
      h <- hist(dat$Hp, breaks=100, plot=F)
      clr <- ifelse(h$breaks > 0-4, "grey", "red")[-length(h$breaks)]
      h <- hist(dat$Hp, breaks=100, freq=T, xlim=xlim, ylim=c(0,max(h$counts)*1.2), col=clr, ylab="", xlab="", main="", plot=T)
      title(main=maintxt, adj=0, outer=F)
      box()
      text(h$breaks[length(h$breaks)*0.25], max(h$counts)*1.1, substitute( mu == A, list(A=round(mean(dat$Hp),3))), adj=c(0,0.5), cex=0.7)
      text(h$breaks[length(h$breaks)*0.75], max(h$counts)*1.1, substitute( sigma == A, list(A=round(sd(dat$Hp),3))), adj=c(0,0.5), cex=0.7)
      mtext(expression(italic(H[P])), side=1, line=2, cex=0.6)
      mtext("n Windows", side=2, line=2, cex=0.6)
    } else {     
      xlim = c(quantile(dat$ZHp,0), quantile(dat$ZHp,1))
      h <- hist(dat$ZHp, breaks=100, plot=F)
      clr <- ifelse(h$breaks > 0-ZLower & h$breaks < ZLower, "grey", "red")[-length(h$breaks)]
      h <- hist(dat$ZHp, breaks=100, freq=T, xlim=xlim, ylim=c(0,max(h$counts)*1.2), col=clr, ylab="", xlab="", main="", plot=T)
      box()
      text(h$breaks[length(h$breaks)*0.25], max(h$counts)*1.1, substitute( mu == A, list(A=round(mean(dat$ZHp),3))), adj=c(0,0.5), cex=0.7)
      text(h$breaks[length(h$breaks)*0.75], max(h$counts)*1.1, substitute( sigma == A, list(A=round(sd(dat$ZHp),3))), adj=c(0,0.5), cex=0.7)
      abline(v=0-ZUpper, lwd=1, lty=3, col=pcolz1)
      abline(v=0-ZLower, lwd=1, lty=3, col=pcolz2)
      abline(v=0+ZUpper, lwd=1, lty=3, col=pcolz1)
      abline(v=0+ZLower, lwd=1, lty=3, col=pcolz2)
      mtext(expression(italic(Z(H[P]))), side=1, line=2, cex=0.6)
      mtext("n Windows", side=2, line=2, cex=0.6)
    }
  }
  dev.off()
}



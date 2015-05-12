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
pdfname <- paste(outpath, "/Results/fst-histograms.png", sep="")
png(pdfname, height=6, width=10, units="in", res=300)
vrow <- 2
vcol <- 4
par(mfrow=c(vrow,vcol), oma = c(1,1,1,1), mar = c(3,3,3,3), cex=0.7, cex.axis=0.7)

for(i in 1:8){
  if(i==1) { win_size <- "1kb"; size <- 1000; maintxt <- "A" }
  if(i==2) { win_size <- "2kb"; size <- 2000; maintxt <- "B" }
  if(i==3) { win_size <- "3kb"; size <- 3000; maintxt <- "C" }
  if(i==4) { win_size <- "4kb"; size <- 4000; maintxt <- "D" }
  if(i==5) { win_size <- "5kb"; size <- 5000; maintxt <- "E" }
  if(i==6) { win_size <- "10kb"; size <- 10000; maintxt <- "F" }
  if(i==7) { win_size <- "20kb"; size <- 20000; maintxt <- "G" }
  if(i==8) { win_size <- "50kb"; size <- 20000; maintxt <- "H" }
  
  fname.fst <- paste("VCFT-fst/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle.windowed.weir.fst", sep="")
  dat.fst <- read.table(fname.fst, header=T, stringsAsFactors=F)
  xlim = c(quantile(dat.fst$N_VARIANTS,0), quantile(dat.fst$N_VARIANTS,1))
  h <- hist(dat.fst$N_VARIANTS, breaks=50, plot=F)
  clr <- ifelse(h$breaks > 10, "grey", "red")[-length(h$breaks)]
  h <- hist(dat.fst$N_VARIANTS, breaks=50, freq=T, xlim=xlim, ylim=c(0,max(h$counts)*1.5), col=clr, ylab="", xlab="", main="", plot=T)
  title(main=maintxt, adj=0, outer=F)
  box()
  text(0, max(h$counts)*1.48, paste("window size: ", win_size, sep=""), adj=c(0,0.5), cex=0.7)
  text(0, max(h$counts)*1.38, paste("number of windows: ", nrow(dat.fst), sep=""), adj=c(0,0.5), cex=0.7)
  text(0, max(h$counts)*1.28, paste("windows with <10 SNPs: ", nrow(dat.fst[which(dat.fst$N_VARIANTS<10),]), sep=""), adj=c(0,0.5), cex=0.7)
  text(0, max(h$counts)*1.18, paste("mean SNPs/window: ", round(mean(dat.fst$N_VARIANTS),2), "+/-", round(sd(dat.fst$N_VARIANTS),2), sep=""), adj=c(0,0.5), cex=0.7)  
  text(0, max(h$counts)*1.08, paste("fraction coverage lost: ", round(nrow(dat.fst[which(dat.fst$N_VARIANTS<10),]) / nrow(dat.fst),2), sep=""), adj=c(0,0.5), cex=0.7)  
  mtext("n SNPs", side=1, line=2, cex=0.6)
  mtext("n Windows", side=2, line=2, cex=0.6)
  
}
dev.off()



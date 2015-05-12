#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(IRanges)
library(stats)
library(pracma)
library(fdrtool)
library(mnormt)

# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/selection")
win_size <- "1kb"
snp_min <- 1
outpath <- "/home/dave/Documents/SeqApiPop"
dir.create(file.path(outpath, paste("Results/", win_size, "/GeneSets", sep="")), recursive=T)

# Clear all data frames
chr <- 0
start <- 0
Val <- 0
dat.hp <- data.frame(chr, start, Val, stringsAsFactors=F)
dat.hp.jfm <- data.frame(chr, start, Val, stringsAsFactors=F)
dat.li.jfm <- data.frame(chr, start, Val, stringsAsFactors=F)
dat.hp.oth <- data.frame(chr, start, Val, stringsAsFactors=F)
dat.li.oth <- data.frame(chr, start, Val, stringsAsFactors=F)
dat.fst <- data.frame(chr, start, Val, stringsAsFactors=False)
dat.daf <- data.frame(chr, start, Val, stringsAsFactors=False)
dat.xpehh <- data.frame(chr, start, Val, stringsAsFactors=False)
dat.pi.oth <- data.frame(chr, start, Val, stringsAsFactors=False)
dat.pi.jfm <- data.frame(chr, start, Val, stringsAsFactors=False)
dat.ihs.oth <- data.frame(chr, start, Val, stringsAsFactors=False)
dat.ihs.jfm <- data.frame(chr, start, Val, stringsAsFactors=False)
dat.alpha.jfm <- data.frame(chr, start, Val, stringsAsFactors=False)
dat.alpha.oth <- data.frame(chr, start, Val, stringsAsFactors=False)




# =========================================================================================================
# Load data
# =========================================================================================================


# Pooled heterozygosity combined data
fname.hp <- paste("Hp/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle.vcf.hp", sep="")
dat.hp <- read.table(fname.hp, header=T, stringsAsFactors=F)
dat.hp <- dat.hp[which(dat.hp$n_snps > snp_min),]
dat.hp <- dat.hp[which(dat.hp$chr<17),]
# Re-calculate Z_Hp (zero mean and unit variance) across all autosomes rather than per chromosome as in file
dat.hp[,11] <- ( dat.hp$Hp - mean(dat.hp$Hp, na.rm=T) ) / sd(dat.hp$Hp, na.rm=T)
names(dat.hp)[1] <- "chr"
names(dat.hp)[2] <- "start"
names(dat.hp)[11] <- "Val"

# Lindley process on Hp for JFM
dat.li <- dat.hp
names(dat.li)[11] <- "ZHp"
dat.li[,12] <- 2*pnorm(-abs(dat.li$ZHp)) # Converts Z score to P value which already has uniform distribution
names(dat.li)[12] <- "UZHp"
#pc <- 0.01 # Set significance threshold for Linldey process at 0.01%
#lizHi_1 <- round(quantile(dat.li$UZHp, p=1-(pc/100), na.rm=T),2)
lizHi_1 <- quantile(dat.li$UZHp, .95)
dat.li[,13] <- -log10(dat.li$UZHp) - lizHi_1 # Create holder field for storing results of Lindley process
# Perform Linldey process on each chromosome independently to avoid overlapping calculations
# Forward direction
for(i in unique(dat.li$chr)){
  scores <- -log10(dat.li[which(dat.li$chr==i),]$UZHp) - lizHi_1
  h=rep(NA,length(scores))
  h[1]=0   
  for (n in 2:length(scores))  h[n]=max(0,h[n-1]+scores[n])
  dat.li[which(dat.li$chr==i),13] <- h
}
names(dat.li)[13] <- "fVal"
# Reverse direction
dat.li[,14] <- dat.li[,13]
for(i in unique(dat.li$chr)){
  scores <- -log10(dat.li[which(dat.li$chr==i),]$UZHp) - lizHi_1
  scores <- rev(scores)
  h=rep(NA,length(scores))
  h[1]=0   
  for (n in 2:length(scores))  h[n]=max(0,h[n-1]+scores[n])
  dat.li[which(dat.li$chr==i),14] <- rev(h)
}
names(dat.li)[14] <- "rVal"
# Mean of forward and reverse
dat.li[,15] <- rowMeans(dat.li[,13:14])
names(dat.li)[15] <- "Val"

# Standardise local score values (may or may not be sensible thing to do)
dat.li[,13] <- (dat.li[,13] - mean(dat.li[,13], na.rm=T) ) /  sd(dat.li[,13], na.rm=T)
dat.li[,14] <- (dat.li[,14] - mean(dat.li[,14], na.rm=T) ) /  sd(dat.li[,14], na.rm=T)
dat.li[,15] <- (dat.li[,15] - mean(dat.li[,15], na.rm=T) ) /  sd(dat.li[,15], na.rm=T)


# Pooled heterozygosity data JFM
fname.hp.jfm <- paste("Hp/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle-JFM.vcf.hp", sep="")
dat.hp.jfm <- read.table(fname.hp.jfm, header=T, stringsAsFactors=F)
dat.hp.jfm <- dat.hp.jfm[which(dat.hp.jfm$n_snps > snp_min),]
dat.hp.jfm <- dat.hp.jfm[which(dat.hp.jfm$chr<17),]
# Re-calculate Z_Hp (zero mean and unit variance) across all autosomes rather than per chromosome as in file
dat.hp.jfm[,11] <- ( dat.hp.jfm$Hp - mean(dat.hp.jfm$Hp, na.rm=T) ) / sd(dat.hp.jfm$Hp, na.rm=T)
names(dat.hp.jfm)[1] <- "chr"
names(dat.hp.jfm)[2] <- "start"
names(dat.hp.jfm)[11] <- "Val"


# Lindley process on Hp for JFM
dat.li.jfm <- dat.hp.jfm
names(dat.li.jfm)[11] <- "ZHp"
dat.li.jfm[,12] <- 2*pnorm(-abs(dat.li.jfm$ZHp)) # Converts Z score to P value which already has uniform distribution
names(dat.li.jfm)[12] <- "UZHp"
#pc <- 0.01 # Set significance threshold for Linldey process at 0.01%
#JFMlizHi_1 <- round(quantile(dat.li.jfm$UZHp, p=1-(pc/100), na.rm=T),2)
JFMlizHi_1 <- quantile(dat.li.jfm$UZHp, .95)
dat.li.jfm[,13] <- -log10(dat.li.jfm$UZHp) - JFMlizHi_1 # Create holder field for storing results of Lindley process
# Perform Linldey process on each chromosome independently to avoid overlapping calculations
# Forward direction
for(i in unique(dat.li.jfm$chr)){
  scores <- -log10(dat.li.jfm[which(dat.li.jfm$chr==i),]$UZHp) - JFMlizHi_1
  h=rep(NA,length(scores))
  h[1]=0   
  for (n in 2:length(scores))  h[n]=max(0,h[n-1]+scores[n])
  dat.li.jfm[which(dat.li.jfm$chr==i),13] <- h
}
names(dat.li.jfm)[13] <- "fVal"
# Reverse direction
dat.li.jfm[,14] <- dat.li.jfm[,13]
for(i in unique(dat.li.jfm$chr)){
  scores <- -log10(dat.li.jfm[which(dat.li.jfm$chr==i),]$UZHp) - JFMlizHi_1
  scores <- rev(scores)
  h=rep(NA,length(scores))
  h[1]=0   
  for (n in 2:length(scores))  h[n]=max(0,h[n-1]+scores[n])
  dat.li.jfm[which(dat.li.jfm$chr==i),14] <- rev(h)
}
names(dat.li.jfm)[14] <- "rVal"
# Mean of forward and reverse
dat.li.jfm[,15] <- rowMeans(dat.li.jfm[,13:14])
names(dat.li.jfm)[15] <- "Val"

# Standardise local score values (may or may not be sensible thing to do)
dat.li.jfm[,13] <- (dat.li.jfm[,13] - mean(dat.li.jfm[,13], na.rm=T) ) /  sd(dat.li.jfm[,13], na.rm=T)
dat.li.jfm[,14] <- (dat.li.jfm[,14] - mean(dat.li.jfm[,14], na.rm=T) ) /  sd(dat.li.jfm[,14], na.rm=T)
dat.li.jfm[,15] <- (dat.li.jfm[,15] - mean(dat.li.jfm[,15], na.rm=T) ) /  sd(dat.li.jfm[,15], na.rm=T)



# Pooled heterozygosity data OTH
fname.hp.oth <- paste("Hp/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle-OTH.vcf.hp", sep="")
dat.hp.oth <- read.table(fname.hp.oth, header=T, stringsAsFactors=F)
dat.hp.oth <- dat.hp.oth[which(dat.hp.oth$n_snps > snp_min),]
dat.hp.oth <- dat.hp.oth[which(dat.hp.oth$chr<17),]
# Re-calculate Z_Hp across all autosomes rather than per chromosome as in file
dat.hp.oth[,11] <- (dat.hp.oth$Hp - mean(dat.hp.oth$Hp, na.rm=T) ) /  sd(dat.hp.oth$Hp, na.rm=T)
names(dat.hp.oth)[1] <- "chr"
names(dat.hp.oth)[2] <- "start"
names(dat.hp.oth)[11] <- "Val"

# Lindley process on Hp for OTH
dat.li.oth <- dat.hp.oth
names(dat.li.oth)[11] <- "ZHp"
dat.li.oth[,12] <- 2*pnorm(-abs(dat.li.oth$ZHp)) # Converts Z score to P value which already has uniform distribution
names(dat.li.oth)[12] <- "UZHp"
#pc <- 0.01 # Set significance threshold for Linldey process at 0.01%
#OTHlizHi_1 <- round(quantile(dat.li.oth$UZHp, p=1-(pc/100), na.rm=T),2)
OTHlizHi_1 <- quantile(dat.li.oth$UZHp, .95)
dat.li.oth[,13] <- -log10(dat.li.oth$UZHp) - OTHlizHi_1 # Create holder field for storing results of Lindley process
# Perform Linldey process on each chromosome independently to avoid overlapping calculations
# Forward direction
for(i in unique(dat.li.oth$chr)){
  scores <- -log10(dat.li.oth[which(dat.li.oth$chr==i),]$UZHp) - OTHlizHi_1
  h=rep(NA,length(scores))
  h[1]=0   
  for (n in 2:length(scores))  h[n]=max(0,h[n-1]+scores[n])
  dat.li.oth[which(dat.li.oth$chr==i),13] <- h
}
names(dat.li.oth)[13] <- "fVal"
# Reverse direction
dat.li.oth[,14] <- dat.li.oth[,13]
for(i in unique(dat.li.oth$chr)){
  scores <- -log10(dat.li.oth[which(dat.li.oth$chr==i),]$UZHp) - OTHlizHi_1
  scores <- rev(scores)
  h=rep(NA,length(scores))
  h[1]=0   
  for (n in 2:length(scores))  h[n]=max(0,h[n-1]+scores[n])
  dat.li.oth[which(dat.li.oth$chr==i),14] <- rev(h)
}
names(dat.li.oth)[14] <- "rVal"
# Mean of forward and reverse
dat.li.oth[,15] <- rowMeans(dat.li.oth[,13:14])
names(dat.li.oth)[15] <- "Val"

# Standardise local score values (may or may not be sensible thing to do)
dat.li.oth[,13] <- (dat.li.oth[,13] - mean(dat.li.oth[,13], na.rm=T) ) /  sd(dat.li.oth[,13], na.rm=T)
dat.li.oth[,14] <- (dat.li.oth[,14] - mean(dat.li.oth[,14], na.rm=T) ) /  sd(dat.li.oth[,14], na.rm=T)
dat.li.oth[,15] <- (dat.li.oth[,15] - mean(dat.li.oth[,15], na.rm=T) ) /  sd(dat.li.oth[,15], na.rm=T)



# VCFT-Fst data
fname.fst <- paste("VCFT-fst/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle.windowed.weir.fst", sep="")
dat.fst <- read.table(fname.fst, header=T, stringsAsFactors=F)
dat.fst <- dat.fst[which(dat.fst$N_VARIANTS > snp_min),]
names(dat.fst)[1] <- "chr"
names(dat.fst)[2] <- "start"
names(dat.fst)[5] <- "Val"


# Comparing difference in windows
#tmp <- dat[,1:3]
#tmp[,4] <- unlist(lapply(seq(1, nrow(tmp)), function(x) paste(tmp[x,1:3],collapse=":")))
#tmp2 <- dat.hp.oth[,1:3]
#tmp2[,4] <- unlist(lapply(seq(1, nrow(tmp2)), function(x) paste(tmp2[x,1:3],collapse=":")))

#DAF
#fname.daf <- "DAF/S2_JFM-OTH.daf.win"
#dat.daf <- read.table(fname.daf, header=T, stringsAsFactors=F)
#names(dat.daf)<- c("chr", "start", "stop", "uDAF", "Val")

# XP-EHH data
#fname.xpehh <- "/home/dave/Copy/HoneyBee/Analyses/selection/xpehh/JFM-OTH.xpehh"
#dat.xpehh <- read.table(fname.xpehh, header=F, stringsAsFactors=F)
#dat.xpehh[,10] <- (dat.xpehh[,9]-mean(dat.xpehh[,9]))/sd(dat.xpehh[,9])
#names(dat.xpehh) <- c("chr", "locus", "start", "gpos", "p1", "ihh1", "p2", "ihh2", "xpehh", "Val")

# VCF-Nucleotide Diversity JFM
fname.pi.jfm <- paste("Nuc-Div/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle-JFM.windowed.pi", sep="")
dat.pi.jfm <- read.table(fname.pi.jfm, header=T, stringsAsFactors=F)
dat.pi.jfm <- dat.pi.jfm[which(dat.pi.jfm$N_VARIANTS > snp_min),]
names(dat.pi.jfm)[1] <- "chr"
names(dat.pi.jfm)[2] <- "start"
names(dat.pi.jfm)[5] <- "Val"

# VCF-Nucleotide Diversity OTH
fname.pi.oth <- paste("Nuc-Div/", win_size, "/JFMOTH_preImp_chrs16_qc_beagle-OTH.windowed.pi", sep="")
dat.pi.oth <- read.table(fname.pi.oth, header=T, stringsAsFactors=F)
dat.pi.oth <- dat.pi.oth[which(dat.pi.oth$N_VARIANTS > snp_min),]
names(dat.pi.oth)[1] <- "chr"
names(dat.pi.oth)[2] <- "start"
names(dat.pi.oth)[5] <- "Val"




# iHS data JFM
#fname.ihs.jfm <- "/home/dave/Copy/HoneyBee/Analyses/selection/iHS/JFM.ihs"
#dat.ihs.jfm <- read.table(fname.ihs.jfm, header=F, stringsAsFactors=F)
#names(dat.ihs.jfm) <- c("chr", "locus", "start", "1-freq", "ihh1", "ihh0", "uniHS", "stdiHS", "Val")
# Two-sided P-value from Z test (neutrality test)
#dat.ihs.jfm[,9] <- -log10(1-2*abs(pnorm(dat.ihs.jfm$stdiHS)-0.5))
# Or Z score
#dat.ihs.jfm[,9] <- dat.ihs.jfm[,8]

# iHS data OTH
#fname.ihs.oth <- "/home/dave/Copy/HoneyBee/Analyses/selection/iHS/OTH.ihs"
#dat.ihs.oth <- read.table(fname.ihs.oth, header=F, stringsAsFactors=F)
#names(dat.ihs.oth) <- c("chr", "locus", "start", "1-freq", "ihh1", "ihh0", "uniHS", "stdiHS", "Val")
# Two-sided P-value from Z test (neutrality test)
#dat.ihs.oth[,9] <- -log10(1-2*abs(pnorm(dat.ihs.oth$stdiHS)-0.5))
# Or Z score
#dat.ihs.oth[,9] <- dat.ihs.oth[,8]

# alpha data JFM
#fname.alpha.jfm <- "/home/dave/Copy/HoneyBee/Analyses/selection/CDS/JFM.alpha"
#dat.alpha.jfm <- read.table(fname.alpha.jfm, header=T, stringsAsFactors=F, sep="\t")
#dat.alpha.jfm[,15] <- -log10(dat.alpha.jfm[,12])
#names(dat.alpha.jfm)[2] <- "chr"
#names(dat.alpha.jfm)[3] <- "start"
#names(dat.alpha.jfm)[15] <- "Val"

# alpha data OTH
#fname.alpha.oth <- "/home/dave/Copy/HoneyBee/Analyses/selection/CDS/OTH.alpha"
#dat.alpha.oth <- read.table(fname.alpha.oth, header=T, stringsAsFactors=F, sep="\t")
#dat.alpha.oth[,15] <- -log10(dat.alpha.oth[,12])
#names(dat.alpha.oth)[2] <- "chr"
#names(dat.alpha.oth)[3] <- "start"
#names(dat.alpha.oth)[15] <- "Val"


# =========================================================================================================



for(plots in c(1,2,3,5,6,7,8,9,10)){
  
  # Data to plot
  #plots <- 1
  if(plots==1) { dat <- dat.fst; pdfname <- paste("Results/", win_size, "/", win_size, "-Fst.png", sep="") }
  if(plots==2) { dat <- dat.pi.jfm; pdfname <- paste("Results/", win_size, "/", win_size, "-Pi-HN.png", sep="") }
  if(plots==3) { dat <- dat.pi.oth; pdfname <- paste("Results/", win_size, "/", win_size, "-Pi-RJ.png", sep="") }
  if(plots==4) { dat <- dat.daf; pdfname <- paste("Results/", win_size, "/", win_size, "-DAF.png", sep="") }
  if(plots==5) { dat <- dat.hp.jfm; pdfname <- paste("Results/", win_size, "/", win_size, "-Hp-HN.png", sep="") }  
  if(plots==6) { dat <- dat.li.jfm; pdfname <- paste("Results/", win_size, "/", win_size, "-Li-HN.png", sep="") }
  if(plots==7) { dat <- dat.hp.oth; pdfname <- paste("Results/", win_size, "/", win_size, "-Hp-RJ.png", sep="") }  
  if(plots==8) { dat <- dat.li.oth; pdfname <- paste("Results/", win_size, "/", win_size, "-Li-RJ.png", sep="") }
  if(plots==9) { dat <- dat.hp; pdfname <- paste("Results/", win_size, "/", win_size, "-Hp-HNRJ.png", sep="") }
  if(plots==10) { dat <- dat.li; pdfname <- paste("Results/", win_size, "/", win_size, "-Li-HNRJ.png", sep="") }
  if(plots==11) { dat <- dat.ihs.jfm; pdfname <- paste("Results/", win_size, "/", win_size, "-iHS-HN.png", sep="") }  
  if(plots==12) { dat <- dat.ihs.oth; pdfname <- paste("Results/", win_size, "/", win_size, "-iHS-RJ.png", sep="") }  
  if(plots==13) { dat <- dat.xpehh; pdfname <- paste("Results/", win_size, "/", win_size, "-XPEHH-HNRJ.png", sep="") }
  
  # File parameters
  setwd(outpath)
  png(pdfname, height=100, width=175, units="mm", res=300)
  
  # Plot parameters
  bkgr="#FFFFFF"
  frgr="#000000"
  pt.size <- 0.2
  pcolz1=rgb(0.0,0.0,1.0)
  pcolz2=rgb(0.8,0.8,1.0)
  vrow <- 1
  vcol <- 1
  par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, cex=0.3, cex.sub=0.7, cex.main=0.7, cex.axis=0.7, cex.lab=0.7, lwd=1,
      mfrow=c(vrow,vcol), oma = c(5.5,1.5,1.5,1), mar = c(1.5,3,1,1))
  
  # Significance thresholds
  pc <- 0.1 
  zHi_1 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)
  pc <- 0.25
  zHi_2 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)      
  # Adjust thresholds for Z-based results
  if(plots > 4) {
    zHi_1 <- 4
    zHi_2 <- 3.719   # Z=3.719 =~ Q=0.0001 =~ 1 in 9999
  }
  
  # Sort data
  tmp <- subset(dat[order(dat$chr, dat$start), ], )
  
  # Set plot dimension limits
  ymin <- min(tmp$Val, na.rm=T) + (min(tmp$Val, na.rm=T) * 0.25)
  if(plots==1) ymin <- 0
  ymax <- max(tmp$Val, na.rm=T) + (max(tmp$Val, na.rm=T) * 0.25)
  if(plots==6 | plots==8 | plots==10) ymax <- max(tmp$fVal, tmp$rVal, na.rm=T) + (max(tmp$fVal, tmp$rVal, na.rm=T) * 0.25)
  
  # Plot
  plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
       xlim=c(0, as.numeric(row.names(tmp[nrow(tmp),]))), ylim = c(ymin, ymax), xlab = NA, 
       ylab = "", las = 1, pch = 19)
  
  # Annotations
  if(plots==1) {  title(ylab = expression(italic(F[ST])), outer=F, line= 2)
              title (main=expression(bold(paste(italic(F[ST]), " Weir & Cockerham 1984", sep=""))), line = 0.6) }
  if(plots==2) { title(ylab = expression(italic(pi)), outer=F, line= 2.5)
             title (main="Nucleotide diversity (HN)", line = 0.6) }
  if(plots==3) { title(ylab = expression(italic(pi)), outer=F, line= 2.5)
             title (main="Nucleotide diversity (RJ)", line = 0.6) }
  if(plots==4) {  title(ylab = "", outer=F, line= 2)
              title (main="Delta allele frequency", line = 0.6) }
  if(plots==5) {  title(ylab = expression(italic(Z(H[P]))), outer=F, line= 2)
              title (main = expression(bold(paste("Pooled heterozygosity (HN)", sep=""))), line = 0.6) }
  if(plots==6) {  title(ylab = expression(italic(Z(H[L]))), outer=F, line= 2)
              title (main= expression(bold(paste(italic(Z(H[P])), " local score (HN)", sep=""))), line = 0.6) }
  if(plots==7) {  title(ylab = expression(italic(Z(H[P]))), outer=F, line= 2)
              title (main = expression(bold(paste("Pooled heterozygosity (RJ)", sep=""))), line = 0.6) }
  if(plots==8) {  title(ylab = expression(italic(Z(H[L]))), outer=F, line= 2)
              title (main= expression(bold(paste(italic(Z(H[P])), " local score (RJ)", sep=""))), line = 0.6) }
  if(plots==9) { title(ylab = expression(italic(Z(H[P]))), outer=F, line= 2)
             title (main = expression(bold(paste("Pooled heterozygosity (HN + RJ)", sep=""))), line = 0.6) }
  if(plots==10) { title(ylab = expression(italic(Z(H[L]))), outer=F, line= 2)
              title (main= expression(bold(paste(italic(Z(H[P])), " local score (HN + RJ)", sep=""))), line = 0.6) }
  if(plots==11) {  title(ylab = "Z(iHS)", outer=F, line= 2)
               title (main="Integrated Haplotype Score (HN)", line = 0.6) }
  if(plots==12) {  title(ylab = "Z(iHS)", outer=F, line= 2)
               title (main="Integrated Haplotype Score (RJ)", line = 0.6) }
  if(plots==13) {  title(ylab = "Z(XP-EHH)", outer=F, line= 2)
               title(main="XP-EHH (HN | RJ)", line = 0.6) }
  
  
  # All points per chromosome
  nchr <- 17
  col <- rep(rgb(0.8,0.8,0.2), nchr)
  col[seq(1,nchr,2)] <- rgb(0.4,0.4,0.4)
  icol <- 1    
  for(n in seq(1,16)){  
    # Local score data forward and reverse lindley process lines
    if(plots==6 | plots==8 | plots==10){
      with(tmp[which(tmp$chr==n),], points(x=as.numeric(row.names(tmp[which(tmp$chr==n),])), y=rVal, col=rgb(1,0.5,0.5,0.5,0.5), pch = 20, cex=pt.size))
      with(tmp[which(tmp$chr==n),], points(x=as.numeric(row.names(tmp[which(tmp$chr==n),])), y=fVal, col=rgb(0.5,1,0.5,0.5,0.5), pch = 20, cex=pt.size))
    }  else {
      with(tmp[which(tmp$chr==n),], points(x=as.numeric(row.names(tmp[which(tmp$chr==n),])), y=Val, col = col[n], pch = 20, cex=pt.size))
    }
    abline(v=max(as.numeric(row.names(tmp[which(tmp$chr==n),]))), col=rgb(0.1,0.1,0.1), lwd=1, lty=3)
    mtext(n, 1, at=mean(as.numeric(row.names(tmp[which(tmp$chr==n),]))), cex=0.7)
  }
  
  # Significant points
  with(tmp[which(tmp$Val > zHi_2),], points(x = row.names(tmp[which(tmp$Val > zHi_2),]), y=Val, col = pcolz2, pch = 20, cex=pt.size))
  with(tmp[which(tmp$Val > zHi_1),], points(x = row.names(tmp[which(tmp$Val > zHi_1),]), y=Val, col = pcolz1, pch = 20, cex=pt.size))
  with(tmp[which(tmp$Val < 0-zHi_2),], points(x = row.names(tmp[which(tmp$Val < 0-zHi_2),]), y=Val, col = pcolz2, pch = 20, cex=pt.size))
  with(tmp[which(tmp$Val < 0-zHi_1),], points(x = row.names(tmp[which(tmp$Val < 0-zHi_1),]), y=Val, col = pcolz1, pch = 20, cex=pt.size))
  
  # Significant thresholds
  abline(h=zHi_2, col=pcolz2)
  abline(h=zHi_1, col=pcolz1)
  abline(h=0-zHi_2, col=pcolz2)
  abline(h=0-zHi_1, col=pcolz1)
  
  # Mean of data
  abline(h=round(mean(dat$Val, na.rm=T),3), col=rgb(0.3,0.3,0.3), lty=2, lwd=1)    
  
  
  box()
  
  dev.off()

}


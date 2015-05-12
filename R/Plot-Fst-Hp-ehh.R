library(stats)
library(pracma)

# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/selection")

# Pooled heterozygosity data JFM
fname.hp.jfm <- "Hp/S2-07112014-JFM.vcf.hp" 
dat.hp.jfm <- read.table(fname.hp.jfm, header=T, stringsAsFactors=F)
dat.hp.jfm <- dat.hp.jfm[which(dat.hp.jfm$n_snps>9),]
dat.hp.jfm <- dat.hp.jfm[which(dat.hp.jfm$chr<17),]
# Re-calculate Z_Hp across all autosomes rather than per chromosome as in file
dat.hp.jfm[,11] <- ( dat.hp.jfm$Hp - mean(dat.hp.jfm$Hp, na.rm=T) ) / sd(dat.hp.jfm$Hp, na.rm=T)
names(dat.hp.jfm)[1] <- "chr"
names(dat.hp.jfm)[2] <- "start"
names(dat.hp.jfm)[11] <- "Val"

# Lindley process on Hp for JFM
dat.li.jfm <- dat.hp.jfm
names(dat.li.jfm)[11] <- "ZHp"
#dat.li.jfm[,12] <- erfc(-dat.li.jfm$ZHp/sqrt(2))/2 # Converts normal distribution (Z) to uniform distribution (U)
dat.li.jfm[,12] <- 2*pnorm(-abs(dat.li.jfm$ZHp)) # Converts Z score to P value which already has uniform distribution
names(dat.li.jfm)[12] <- "UZHp"
pc <- 0.01 # Set significance threshold for Linldey process at 0.01%
JFMlizHi_1 <- round(quantile(dat.li.jfm$UZHp, p=1-(pc/100), na.rm=T),2)
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



# Pooled heterozygosity data OTH
fname.hp.oth <- "Hp/S2-07112014-OTH.vcf.hp"
dat.hp.oth <- read.table(fname.hp.oth, header=T, stringsAsFactors=F)
dat.hp.oth <- dat.hp.oth[which(dat.hp.oth$n_snps>9),]
dat.hp.oth <- dat.hp.oth[which(dat.hp.oth$chr<17),]
# Re-calculate Z_Hp across all autosomes rather than per chromosome as in file
dat.hp.oth[,11] <- (dat.hp.oth$Hp - mean(dat.hp.oth$Hp, na.rm=T) ) /  sd(dat.hp.oth$Hp, na.rm=T)
names(dat.hp.oth)[1] <- "chr"
names(dat.hp.oth)[2] <- "start"
names(dat.hp.oth)[11] <- "Val"

# Lindley process on Hp for OTH
dat.li.oth <- dat.hp.oth
names(dat.li.oth)[11] <- "ZHp"
#dat.li.oth[,12] <- erfc(-dat.li.oth$ZHp/sqrt(2))/2 # Converts normal distribution (Z) to uniform distribution (U)
dat.li.oth[,12] <- 2*pnorm(-abs(dat.li.oth$ZHp)) # Converts Z score to P value which already has uniform distribution
names(dat.li.oth)[12] <- "UZHp"
pc <- 0.01 # Set significance threshold for Linldey process at 0.01%
OTHlizHi_1 <- round(quantile(dat.li.oth$UZHp, p=1-(pc/100), na.rm=T),2)
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



# VCFT-Fst data
fname.fst <- "VCFT-fst/S2-07112014_hapJFMOTH.windowed.weir.fst"
dat.fst <- read.table(fname.fst, header=T, stringsAsFactors=F)
dat.fst <- dat.fst[which(dat.fst$N_VARIANTS>9),]
names(dat.fst)[1] <- "chr"
names(dat.fst)[2] <- "start"
names(dat.fst)[5] <- "Val"

# VCF-Nucleotide Diversity JFM
#fname.pi.jfm <- "Nuc-Div/JFM.windowed.pi"
#dat.pi.jfm <- read.table(fname.pi.jfm, header=T, stringsAsFactors=F)
#dat.pi.jfm <- dat.pi.jfm[which(dat.pi.jfm$N_VARIANTS>9),]
#names(dat.pi.jfm)[1] <- "chr"
#names(dat.pi.jfm)[2] <- "start"
#names(dat.pi.jfm)[5] <- "Val"

# VCF-Nucleotide Diversity OTH
#fname.pi.oth <- "Nuc-Div/OTH.windowed.pi"
#dat.pi.oth <- read.table(fname.pi.oth, header=T, stringsAsFactors=F)
#dat.pi.oth <- dat.pi.oth[which(dat.pi.oth$N_VARIANTS>9),]
#names(dat.pi.oth)[1] <- "chr"
#names(dat.pi.oth)[2] <- "start"
#names(dat.pi.oth)[5] <- "Val"

#DAF
#fname.daf <- "DAF/S2_JFM-OTH.daf.win"
#dat.daf <- read.table(fname.daf, header=T, stringsAsFactors=F)
#names(dat.daf)<- c("chr", "start", "stop", "uDAF", "Val")

# XP-EHH data
#fname.xpehh <- "/home/dave/Copy/HoneyBee/Analyses/selection/ehh/D-JFM-OTH-all.xpehh.out"
#dat.xpehh <- read.table(fname.xpehh, header=T, stringsAsFactors=F)
#dat.xpehh[,9] <- (dat.xpehh[,8]-mean(dat.xpehh[,8]))/sd(dat.xpehh[,8])
#names(dat.xpehh)[2] <- "start"
#names(dat.xpehh)[9] <- "Val"


# Add chromosome and re-label vields
plot_x_label <- ""
n <- 1
n_plots <- 5

for(n in seq(1,16))
{
  
  chromosome <- n # Set to 0 plots chromosomes 1-16
  savetopdf <- 1
  
  
  # Plot parameters
  if(savetopdf==1){
  #  pdfname <- paste(chromosome, ".pdf", sep="")
  #  pdf(pdfname, height=9, width=6)
    pdfname <- paste(chromosome, ".png", sep="")
    png(pdfname, height=9, width=6, units="in", res=300)
  }
  
  bkgr="#FFFFFF"
  frgr="#000000"
  pcolz1=rgb(0.0,0.0,1.0)
  pcolz2=rgb(0.8,0.8,1.0)
  par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
      cex=0.3, cex.sub=0.7, cex.main=0.7, cex.axis=0.7, cex.lab=1, lwd=1)
  vrow <- n_plots
  vcol <- 1
  par(mfrow=c(vrow,vcol), oma = c(5.5,2.5,1,1), mar = c(1.5,0.5,1,0))
  
  
  
  # ======================================================================================================
  # PLOT
  # ======================================================================================================
  
  for(i in seq(1,n_plots)){

    # 3 and 5 MUST be left for dat.li.* due to hard-coded plot differences
  if(i==1) dat <- dat.fst
  if(i==2) dat <- dat.hp.jfm
  if(i==3) dat <- dat.li.jfm
  if(i==4) dat <- dat.hp.oth
  if(i==5) dat <- dat.li.oth
  if(i==6) dat <- dat.xpehh
  if(i==7) dat <- dat.daf  
  
  # Quantile
  pc <- 0.01 # 0.01%
  zHi_1 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)
  pc <- 0.1 # 0.1%
  zHi_2 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)

  # Sort data and get max and min Z_Hp values
  tmp <- dat[which(dat$chr==chromosome),]
  tmp <- subset(tmp[order(tmp$chr, tmp$start), ], )
  ymax <- (max(tmp$Val, na.rm=T)) + (0.2*max(tmp$Val, na.rm=T))
  ymin <- (min(tmp$Val, na.rm=T)) - (0.2*max(tmp$Val, na.rm=T))
  if(i==3 | i==5){
    ymax <- (max(tmp$fVal, na.rm=T)) + (0.2*max(tmp$fVal, na.rm=T))
    ymin <- (min(tmp$fVal, na.rm=T)) - (0.2*max(tmp$fVal, na.rm=T))
  }
  
  xmin <- 0
  xmax <- max(dat.fst[which(dat.fst$chr==chromosome),]$start, 
              dat.hp.jfm[which(dat.hp.jfm$chr==chromosome),]$start,
              dat.hp.oth[which(dat.hp.oth$chr==chromosome),]$start,
              dat.li.jfm[which(dat.li.jfm$chr==chromosome),]$start,
              dat.li.oth[which(dat.li.oth$chr==chromosome),]$start)
  #dat.pi.jfm[which(dat.pi.jfm$chr==chromosome),]$start,
  #dat.pi.oth[which(dat.pi.oth$chr==chromosome),]$start,
#              dat.daf[which(dat.daf$chr==chromosome),]$start,
#              dat.xpehh[which(dat.xpehh$chr==chromosome),]$start)
  
  # Plot
  plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
       xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = NA, 
       ylab = plot_x_label, las = 1, pch = 19)

  # Line for fVal and rVal from Lindley process
  if(i==3 || i==5){
    with(tmp, points(x=start, y=tmp$fVal, type="l", lwd=1, col=rgb(1,0.5,0.5)))
    with(tmp, points(x=start, y=tmp$rVal, type="l", lwd=1, col=rgb(0.5,1,0.5)))
  }

  # Plot points for all types except for Lindley process
  nchr <- 17
  col <- rep(rgb(0.6,0.6,0.6), nchr)
  col[seq(1,nchr,2)] <- rgb(0.8,0.8,0.8)
  icol <- 1    

  # Significance thresholds and mean indicator
  if(i!=3 & i!=5) {
    with(tmp, points(x=start, y=Val, col = col[icol], pch = 20))
    with(tmp[which(tmp$Val > zHi_2),], points(x=start, y=Val, col = pcolz2, pch = 20))
    with(tmp[which(tmp$Val > zHi_1),], points(x=start, y=Val, col = pcolz1, pch = 20))
    with(tmp[which(tmp$Val < 0-zHi_2),], points(x=start, y=Val, col = pcolz2, pch = 20))
    with(tmp[which(tmp$Val < 0-zHi_1),], points(x=start, y=Val, col = pcolz1, pch = 20))
    abline(h=zHi_1, col=pcolz1)
    abline(h=zHi_2, col=pcolz2)
    abline(h=0-zHi_1, col=pcolz1)
    abline(h=0-zHi_2, col=pcolz2)
    abline(h=round(mean(dat$Val, na.rm=T),3), col=rgb(0.3,0.3,0.3), lty=2, lwd=1)    
  } else {    
    with(tmp, points(x=start, y=Val,  col = col[icol], type="l", lwd=2))
    pc <- 5
    if(i==3) zHi_1 <- round(quantile(dat.li.jfm$Val, p=1-(pc/100), na.rm=T),2)
    if(i==5) zHi_1 <- round(quantile(dat.li.oth$Val, p=1-(pc/100), na.rm=T),2)    
    abline(h=zHi_1, col=pcolz1, lwd=1)
    tmp[which(tmp$Val< zHi_1),]$Val <- NA
    with(tmp, points(x=start, y=Val, type="l", lwd = 3, col = pcolz1))
  }

  
  if(i==1) title (main="Fst [Weir & Cockerham 1984]", line = 0.4)
  if(i==2) title (main="Pooled heterozygosity (JFM)", line = 0.4)
  if(i==3) title (main="Lindley equation on Hp (JFM)", line = 0.4)
  if(i==4) title (main="Pooled heterozygosity (OTH)", line = 0.4)
  if(i==5) title (main="Lindley equation on Hp (OTH)", line = 0.4)
  if(i==6) title (main="Nucleotide diversity (JFM)", line = 0.4)
  if(i==7) title (main="Nucleotide diversity (OTH)", line = 0.4)
  if(i==8) title(main="XP-EHH (JFM | OTH)", line = 0.4)
  if(i==9) title (main="Delta allele frequency", line = 0.4)
  if(i==vrow) {
    title(xlab = "Position", outer=T)
    axis(1, col=frgr, at=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))],
         labels=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))]) 
  }
  
  
  
    # Known loci (not apparently published, are from a Chinese group's presentation)
    ltype=2
    lwidth=1
    # ROYAL JELLY
    # Dop2      LG15:6784429..6844358
    # hex7  	  LG2:12690182..12694059
    # SsRbeta		LG10:13320..15501
    hex7 <- tmp[which(tmp$chr==2)[1],]$start + 12690182
    SsRbeta <- tmp[which(tmp$chr==10)[1],]$start + 13320
    Dop2 <- tmp[which(tmp$chr==15)[1],]$start + 6784429
#    abline(v=hex7, col=rgb(1.0,0.5,1.0), lty=ltype, lwd=lwidth)
#    text(x=hex7, y=(ymin*0.9), labels="hex7", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))
#    abline(v=SsRbeta, col=rgb(1.0,0.5,1.0), lty=ltype, lwd=lwidth)
#    text(x=SsRbeta, y=(ymin*0.9), labels="SsRbeta", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))
#    abline(v=Dop2, col=rgb(1.0,0.5,1.0), lty=ltype, lwd=lwidth)
#    text(x=Dop2, y=(ymin*0.9), labels="Dop2", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))
    
    # MRJP genes (all between MRJPL and MRJPR)
    MRJPL <-   tmp[which(tmp$chr==11)[1],]$start + 2555000
    abline(v=MRJPL, col=rgb(0.5,0.5,1.0), lty=ltype, lwd=lwidth)
    MRJPR <-   tmp[which(tmp$chr==11)[1],]$start + 2637000
    abline(v=MRJPR, col=rgb(0.5,0.5,1.0), lty=ltype, lwd=lwidth)
        
    # FORAGING
    # For       LG13:8976929..9056381
    For <-  tmp[which(tmp$chr==13)[1],]$start + 8976929
    abline(v=For, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
    text(x=For, y=(ymin*0.9), labels="For", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))
    # QTL from Page et al. 2012 - Complex pleiotropy...
    # stsD8 (pln1) QTL  LG13:3.4Mb
    stsD8 <-  tmp[which(tmp$chr==13)[1],]$start + 3400000
    abline(v=stsD8, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
    text(x=stsD8, y=(ymin*0.9), labels="pln1", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))
    # E15 (pln1) QTL  LG13:7.1Mb
    E15 <-  tmp[which(tmp$chr==13)[1],]$start + 7100000
    abline(v=E15, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
    text(x=E15, y=(ymin*0.9), labels="pln1", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))
    # W5 (pln2) QTL  LG1:16.9Mb
    W5 <-  tmp[which(tmp$chr==1)[1],]$start + 16900000
    abline(v=W5, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
    text(x=W5, y=(ymin*0.9), labels="pln2", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))
    # F8 (pln2) QTL  LG1:19.7Mb
    F8 <-  tmp[which(tmp$chr==1)[1],]$start + 19700000
    abline(v=F8, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
    text(x=F8, y=(ymin*0.9), labels="pln2", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))
    # stsQ4 (pln3) QTL  LG1:9.2Mb
    stsQ4 <-  tmp[which(tmp$chr==1)[1],]$start + 9200000
    abline(v=stsQ4, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
    text(x=stsQ4, y=(ymin*0.9), labels="pln3", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))
  
    # Compelementary sex-determining locus (csd) LG3:11135135-11146488
    stsQ4 <-  tmp[which(tmp$chr==3)[1],]$start + 11140812
    abline(v=stsQ4, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
    text(x=stsQ4, y=(ymin*0.9), labels="CSD", pos=4, offset=0.5,cex=0.6,col=rgb(0.0,0.0,1.0))

    
  # Box-it
  box()
  
  }
  
  
  
  if(savetopdf==1) dev.off()

}





# write.table(dat.fst[which(dat.fst$Val > zHi_1),c(1:3,5),], "VCFT-fst/fst.bed", sep="\t", row.names=F, col.names=F, append=F, quote=F)
# write.table(dat.hp.jfm[which(dat.hp.jfm$Val < 0-3.34),c(1:3,11),], "Hp/JFM-neg.bed", sep="\t", row.names=F, col.names=F, append=F, quote=F)
# write.table(dat.hp.jfm[which(dat.hp.jfm$Val > 3.34),c(1:3,11),], "Hp/JFM-pos.bed", sep="\t", row.names=F, col.names=F, append=F, quote=F)
# write.table(dat.hp.oth[which(dat.hp.oth$Val < 0-2.9),c(1:3,11),], "Hp/OTH-neg.bed", sep="\t", row.names=F, col.names=F, append=F, quote=F)
# write.table(dat.hp.oth[which(dat.hp.oth$Val > 2.9),c(1:3,11),], "Hp/OTH-pos.bed", sep="\t", row.names=F, col.names=F, append=F, quote=F)



# Datla normality plots
#png("data-normality.png", height=6, width=6, units="in", res=300)
#bkgr="#FFFFFF"
#frgr="#000000"
#par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
#    cex=0.3, cex.sub=0.7, cex.main=0.7, cex.axis=0.7, cex.lab=1, lwd=1)
#vrow <- 2
#vcol <- 2
#par(mfrow=c(vrow,vcol), oma = c(1,1,1,1), mar = c(2,2,2,2))
#hist(dat.hp.oth$Val)
#qqnorm(dat.hp.oth$Val, main="Normal Q-Q Plot of dat.hp.oth$Val", plot.it=T, datax=F, pch=20)
#qqline(dat.hp.oth$Val)
#hist(dat.hp.jfm$Val)
#qqnorm(dat.hp.jfm$Val, main="Normal Q-Q Plot of dat.hp.jfm$Val", plot.it=T, datax=F, pch=20)
#qqline(dat.hp.jfm$Val)
#dev.off()



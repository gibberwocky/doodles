library(stats)
library(pracma)

# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/selection")

# VCFT-Fst data
fname.fst <- "VCFT-fst/S2-07112014_hapJFMOTH.windowed.weir.fst"
dat.fst1 <- read.table(fname.fst, header=T, stringsAsFactors=F)
dat.fst1 <- dat.fst1[which(dat.fst1$N_VARIANTS>9),]
names(dat.fst1)[1] <- "chr"
names(dat.fst1)[2] <- "start"
names(dat.fst1)[5] <- "Val"
fname.fst <- "VCFT-fst/S2-07112014_hapJFMOTH_preImp.windowed.weir.fst"
dat.fst2 <- read.table(fname.fst, header=T, stringsAsFactors=F)
dat.fst2 <- dat.fst2[which(dat.fst2$N_VARIANTS>9),]
names(dat.fst2)[1] <- "chr"
names(dat.fst2)[2] <- "start"
names(dat.fst2)[5] <- "Val"


# Pooled heterozygosity data JFM
fname.hp.jfm1 <- "Hp/S2-07112014-OTH.vcf.hp" 
dat.hp.jfm1 <- read.table(fname.hp.jfm1, header=T, stringsAsFactors=F)
dat.hp.jfm1 <- dat.hp.jfm1[which(dat.hp.jfm1$n_snps>9),]
dat.hp.jfm1 <- dat.hp.jfm1[which(dat.hp.jfm1$chr<17),]
# Re-calculate Z_Hp across all autosomes rather than per chromosome as in file
dat.hp.jfm1[,11] <- ( dat.hp.jfm1$Hp - mean(dat.hp.jfm1$Hp, na.rm=T) ) / sd(dat.hp.jfm1$Hp, na.rm=T)
names(dat.hp.jfm1)[1] <- "chr"
names(dat.hp.jfm1)[2] <- "start"
names(dat.hp.jfm1)[11] <- "Val"

fname.hp.jfm2 <- "Hp/OTH_qc-AN.vcf.hp" 
dat.hp.jfm2 <- read.table(fname.hp.jfm2, header=T, stringsAsFactors=F)
dat.hp.jfm2 <- dat.hp.jfm2[which(dat.hp.jfm2$n_snps>9),]
dat.hp.jfm2 <- dat.hp.jfm2[which(dat.hp.jfm2$chr<17),]
# Re-calculate Z_Hp across all autosomes rather than per chromosome as in file
dat.hp.jfm2[,11] <- ( dat.hp.jfm2$Hp - mean(dat.hp.jfm2$Hp, na.rm=T) ) / sd(dat.hp.jfm2$Hp, na.rm=T)
names(dat.hp.jfm2)[1] <- "chr"
names(dat.hp.jfm2)[2] <- "start"
names(dat.hp.jfm2)[11] <- "Val"


# Add chromosome and re-label vields
plot_x_label <- ""
n <- 1
n_plots <- 2

for(n in seq(1,16))
{
  
  chromosome <- n # Set to 0 plots chromosomes 1-16
  savetopdf <- 1
  
  
  # Plot parameters
  if(savetopdf==1){
  #  pdfname <- paste(chromosome, ".pdf", sep="")
  #  pdf(pdfname, height=9, width=6)
    pdfname <- paste(chromosome, "_oth_test.png", sep="")
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
  #if(i==1) dat <- dat.fst1
  #if(i==2) dat <- dat.fst2
  if(i==1) dat <- dat.hp.jfm1
  if(i==2) dat <- dat.hp.jfm2
  
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
  xmax <- max(dat.hp.jfm1[which(dat.hp.jfm1$chr==chromosome),]$start, 
              dat.hp.jfm2[which(dat.hp.jfm2$chr==chromosome),]$start)
  
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

  
  #if(i==1) title (main="Fst [Weir & Cockerham 1984] imputed", line = 0.4)
  #if(i==2) title (main="Fst [Weir & Cockerham 1984] non-imputed", line = 0.4)
  if(i==1) title (main="Pooled heterozygosity (OTH) imputed", line = 0.4)
  if(i==2) title (main="Pooled heterozygosity (OTH) non-imputed", line = 0.4)
  
  
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



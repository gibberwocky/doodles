setwd("/home/dave/Copy/HoneyBee/Analyses/selection")

fname <- "1-Hp.out"
dat <- read.table(fname, header=T, stringsAsFactors=F)
pdfname <- paste(fname, ".pdf", sep="")
pdf(pdfname, height=9, width=12)

# Unique chromosome names
chrs <- as.vector(sort(unique(dat$chr)))

# Get size of chromosomes relative to chromosome 1
chrLens <- NULL
for(i in chrs){
  chrLens[i] <- max(dat[which(dat$chr==i),]$stop)
  if(i>1) chrLens[i] <- round(chrLens[i]/chrLens[1],2)
}
chrLens[1] <- 1

# Manual calculation of layout of chromosomes in multiplot
# chrLens = 7.33
# [1.00] [0.52] [0.44] 
# [0.43] [0.48] [0.62] [0.44]
# [0.45] [0.37] [0.43] [0.49]
# [0.40] [0.34] [0.34] [0.34] [0.24]

# 4 2 2
# 2 2 3 2
# 2 2 2 2
# 2 2 2 2 1

vrow <- 4
vcol <- 9
#par(mfrow=c(vrow,vcol),
#    oma = c(2,2,1,1) ,
#    mar = c(0.5,0.5,0,0) )
layout(matrix(c(1, 1, 1, 1, 2, 2, 3, 3, 0,
                4, 4, 5, 5, 6, 6, 6, 7, 7,
                8, 8, 9, 9, 10, 10, 11, 11, 0,
                12, 12, 13, 13, 14, 14, 15, 15, 16),vrow,vcol,byrow=T),
       widths=c(1,1), heights=c(1,1), T)
par(oma = c(4,4,1,1), mar = c(0.5,0.5,0,0) )
par(cex=0.8, cex.main=1, cex.lab=1.2, cex.axis=0.8, cex.sub=1)
par(bg="#333333", col="white", col.axis="white", col.lab="white", col.main="white", fg="white")

# Plots on far left
col1 <- c(1,4,8,12)

for(i in chrs){
  tmp <- dat[which(dat$chr==i),]
  tmp[1,]$Z_Hp <- max(dat$Z_Hp, na.rm=T)
  tmp[2,]$Z_Hp <- min(dat$Z_Hp, na.rm=T)

  if(i%in%col1) {
    plot(tmp$Z_Hp, type="n", ylab="", xlab="", main="", xaxt="n")
  } else {
    plot(tmp$Z_Hp, type="n", ylab="", xlab="", main="",mgp=c(1,0.2,0), yaxt="n", xaxt="n")
  }
  
  tmp <- dat[which(dat$chr==i),]
  tmp[which(tmp$Z_Hp<3 & tmp$Z_Hp>(0-3)),] <- NA
  points(tmp$Z_Hp, type="p", lwd=1.5, pch=20,
         col="#FFD320", ylab="", xlab="", main="")
  tmp <- dat[which(dat$chr==i),]
  tmp[which(tmp$Z_Hp>3),] <- NA
  tmp[which(tmp$Z_Hp<(0-3)),] <- -NA
  points(tmp$Z_Hp, type="p", lwd=1.5, pch=20,
         col=rgb(0.6,0.6,0.6), ylab="", xlab="", main="")
  
  abline(h=3, col=rgb(0.0,0.75,0.0))
  abline(h=-3, col=rgb(0.0,0.75,0.0))
  title(main=paste(" ", i, sep=""), adj=0, line=-1)
  
}

title (xlab="Chromosomal position", outer=T, line=-8)
title (ylab="Z-Hp", outer=T, line=2)
dev.off()



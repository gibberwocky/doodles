# Have run with several different settings (window size, step size) in an 
# effort to find the optimal sized window. This is often subjective. One
# suggestion I think was to go with the size giving the best fit normal
# distribution for the number of SNPs per window. Currently testing for this.
setwd("/home/dave/Copy/HoneyBee/Analyses/selection/OTH")

#fname <- "1-Hp-10k-25step.out"     # A = 515.4015
#fname <- "1-Hp-20k-25step.out"     # A = 342.925
#fname <- "1-Hp-40k-25step.out"     # A = 184.0689
#fname <- "1-Hp-60k-25step.out"     # A = 115.737
#fname <- "1-Hp-80k-25step.out"     # A = 66.2334
fname <- "2-Hp-100k-25step.out"    # A = 38.8456
#fname <- "1-Hp-120k-25step.out"    # A = 23.5
#fname <- "1-Hp-140k-25step.out"    # A = 15.2142
#fname <- "1-Hp-160k-25step.out"    # A = 10.2639
dat <- read.table(fname, header=T, stringsAsFactors=F)

# Re-calculate Z_Hp across all autosomes rather than per chromosome as in file
dat[which(dat$chr < 17),11] <- (dat$Hp - mean(dat[which(dat$chr < 17),]$Hp, na.rm=T) ) / sd(dat[which(dat$chr < 17),]$Hp, na.rm=T)

# Remove windows with fewer than 10 SNPs
dat <- dat[which(dat$n_SNPs>9),]


# ==============================================================================
# Check for normality (Anderson-Darling normality test)
# ==============================================================================
library(nortest)
ad.test(dat$n_SNPs)
ad.values <- c(515, 342, 184, 115, 66, 38, 23, 15, 10)
names(ad.values) <- c("10k", "20k", "40k", "60k", "80k", "100k", "120k", "140k", "160k")
# Calculate mean and SD to identify where variance in AD values slopes off
# Generate plot and mark the window size where variance becomes less than
# the mean of remaining windows + 1 standard deviation
ad.sd <- NULL
ad.mean <- NULL
for(i in 1:length(ad.values)){
  ad.sd[i] <- sd(ad.values[i:length(ad.values)])
  ad.mean[i] <- mean(ad.values[i:length(ad.values)])
}
ad.sd[which(is.na(ad.sd))] <- 0
ad.var <- ad.values <= ad.mean + (1.5*ad.sd)
ad.var[which(ad.var==T)] <- ad.values[ad.var]
# Plot file name
pdfname <- paste(fname, ".adnorm.pdf", sep="")
pdf(pdfname, height=3, width=9)
# Dark background
#bkgr="#000000"
#frgr="#FFFFFF"
#pcol="#FFD320"
# White background
bkgr="#FFFFFF"
frgr="#000000"
pcol=rgb(0.0,0.7,0.0)
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
    cex=0.5, cex.sub=1, cex.main=1, cex.axis=1, cex.lab=1.5, lwd=1.5)
plot(ad.values, pch=19, xaxt="n", ylab="A-value", xlab="Window size", col=pcol)
axis(1, labels=names(ad.values), at=seq(1:length(ad.values)))
abline(v=which(ad.var>0)[1], col=pcol)
dev.off()

# Normality plot for a single window size
# Plot file name
pdfname <- paste(fname, ".qqnorm.pdf", sep="")
pdf(pdfname, height=4, width=4)
par(bg="#333333", fg="#FFFFFF", col.lab="#FFFFFF", col.axis="#FFFFFF", 
    cex=0.5, cex.sub=1, cex.main=1, cex.axis=1, cex.lab=1.5, lwd=1.5)
qqnorm(dat$n_SNPs, main="", pch=20, col="#FFD320")
dev.off()






# ==============================================================================
# Generate pooled heterozygosity plot
#
# Analyse candidate regions on finder scale:
# JFM (100k 25 step) 2, 6, 8, 11
# OTH (100k 25 step) 
# ==============================================================================

# Note significance threshold (0.01% tails)
pc <- 0.01
zQuant <- round(quantile(dat$Z_Hp, p=1-(pc/100), na.rm=T),2)
chromosome <- 0 # Set to 0 plots chromosomes 1-16

# Sort data and get max and min Z_Hp values
tmp <- dat
if(chromosome > 0) tmp <- dat[which(dat$chr==chromosome),]
tmp <- subset(tmp[order(tmp$chr, tmp$start), ], )
ymax <- ceiling(max(tmp$Z_Hp))
ymin <- ceiling(min(tmp$Z_Hp))

# Index data
tmp$index <- NA
ind <- 0
for (i in unique(tmp$chr)) {
  ind = ind + 1
  tmp[tmp$chr == i, ]$index <- ind
}

# Generate 'ticks'
nchr <- length(unique(tmp$chr))
tmp$pos <- NA
lastbase <- 0
ticks <- NULL
for(i in unique(tmp$index)){
  if (i == 1) {
    tmp[tmp$index == i, ]$pos = tmp[tmp$index == i, ]$start
  }
  else {
    lastbase = lastbase + tail(subset(tmp, index == i - 1)$start, 1)
    tmp[tmp$index == i, ]$pos = tmp[tmp$index == i, ]$start + lastbase
  }
  lastbase <- lastbase + 2500000
  ticks = c(ticks, tmp[tmp$index == i, ]$pos[floor(length(tmp[tmp$index == i, ]$pos)/2) + 1])
}
xlabel <- "Chromosome"
labs <- unique(tmp$chr)
xmax = ceiling(max(tmp$pos) * 1.03)
xmin = floor(max(tmp$pos) * -0.03)

# Plot
pdfname <- paste(fname, ".pdf", sep="")
pdf(pdfname, height=3, width=9)

# White background
bkgr="#FFFFFF"
frgr="#000000"
pcol=rgb(0.0,0.8,0.0)
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
    cex=0.5, cex.sub=1, cex.main=1, cex.axis=1, cex.lab=1.5, lwd=1.5)
plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
     xlim = c(xmin, xmax), ylim = c(ymin-2, ymax+2), xlab = xlabel, 
     ylab = paste("Z-Hp (p < 0.01 ~ |Z-Hp| > ", zQuant, ")", sep=""), las = 1, pch = 20)
axis(1, at=ticks, labels=labs, lty=0)

# Points
col <- rep(rgb(0.2,0.2,0.2), nchr)
col[seq(1,nchr,2)] <- rgb(0.5,0.5,0.5)
tmp[which(tmp$Z_Hp >= (0-zQuant) & tmp$Z_Hp <=zQuant & (tmp$chr %% 2 == 0)),3] <- rgb(0.2,0.2,0.2)
icol <- 1
for (i in unique(tmp$index)) {
  with(tmp[tmp$index == unique(tmp$index)[i], ], points(pos, Z_Hp, col = col[icol], pch = 20))
  icol = icol + 1
}
with(tmp[which(tmp$Z_Hp <= (0-zQuant) | tmp$Z_Hp >= zQuant),], points(pos, Z_Hp, col = pcol, pch = 20))

# Significance thresholds
abline(h=zQuant, col=rgb(1.0,0.1,0.0))
abline(h=0-zQuant, col=rgb(1.0,0.1,0.0))
box()

# Known loci
if(chromosome==0){
  # ROYAL JELLY
  # Dop2      LG15:6784429..6844358
  # hex7  	  LG2:12690182..12694059
  # SsRbeta		LG10:13320..15501
  hex7 <- tmp[which(tmp$chr==2)[1],]$pos + 12690182
  SsRbeta <- tmp[which(tmp$chr==10)[1],]$pos + 13320
  Dop2 <- tmp[which(tmp$chr==15)[1],]$pos + 6784429
  abline(v=hex7, col=rgb(0.5,0.5,1.0), lty=3)
  abline(v=SsRbeta, col=rgb(0.5,0.5,1.0), lty=3)
  abline(v=Dop2, col=rgb(0.5,0.5,1.0), lty=3)
  # FORAGING
  # For       LG13:8976929..9056381
  For <-  tmp[which(tmp$chr==13)[1],]$pos + 8976929
  abline(v=For, col=rgb(1.0,0.5,0.0), lty=3)
  # QTL from Page et al. 2012 - Complex pleiotropy...
  # stsD8 (pln1) QTL  LG13:3.4Mb
  stsD8 <-  tmp[which(tmp$chr==13)[1],]$pos + 3400000
  abline(v=stsD8, col=rgb(1.0,0.5,0.0), lty=3)
  # E15 (pln1) QTL  LG13:7.1Mb
  E15 <-  tmp[which(tmp$chr==13)[1],]$pos + 7100000
  abline(v=E15, col=rgb(1.0,0.5,0.0), lty=3)
  # W5 (pln2) QTL  LG1:16.9Mb
  W5 <-  tmp[which(tmp$chr==1)[1],]$pos + 16900000
  abline(v=W5, col=rgb(1.0,0.5,0.0), lty=3)
  # F8 (pln2) QTL  LG1:19.7Mb
  F8 <-  tmp[which(tmp$chr==1)[1],]$pos + 19700000
  abline(v=F8, col=rgb(1.0,0.5,0.0), lty=3)
  # stsQ4 (pln3) QTL  LG1:9.2Mb
  stsQ4 <-  tmp[which(tmp$chr==1)[1],]$pos + 9200000
  abline(v=stsQ4, col=rgb(1.0,0.5,0.0), lty=3)
}



title(main=fname)

# Write to pdf
dev.off()










# ==============================================================================
# Write candidate regions to file
# ==============================================================================

#Subset to hits
hits <- dat[which(dat$Z_Hp<(0-zQuant) | dat$Z_Hp>(zQuant)),1:3]

# Fix chromosome names1 again
hits[which(hits[,1]=="1"),1] <- "NC_007070.3"
hits[which(hits[,1]=="2"),1] <- "NC_007071.3"
hits[which(hits[,1]=="3"),1] <- "NC_007072.3"
hits[which(hits[,1]=="4"),1] <- "NC_007073.3"
hits[which(hits[,1]=="5"),1] <- "NC_007074.3"
hits[which(hits[,1]=="6"),1] <- "NC_007075.3"
hits[which(hits[,1]=="7"),1] <- "NC_007076.3"
hits[which(hits[,1]=="8"),1] <- "NC_007077.3"
hits[which(hits[,1]=="9"),1] <- "NC_007078.3"
hits[which(hits[,1]=="10"),1] <- "NC_007079.3"
hits[which(hits[,1]=="11"),1] <- "NC_007080.3"
hits[which(hits[,1]=="12"),1] <- "NC_007081.3"
hits[which(hits[,1]=="13"),1] <- "NC_007082.3"
hits[which(hits[,1]=="14"),1] <- "NC_007083.3"
hits[which(hits[,1]=="15"),1] <- "NC_007084.3"
hits[which(hits[,1]=="16"),1] <- "NC_007085.3"
hits[which(hits[,1]=="17"),1] <- "NC_001566.1"

# Write as BED file
hitsfname <- paste(fname, ".bed", sep="")
write.table(hits, hitsfname, sep="\t", row.names=F, col.names=F, quote=F)

# Write as BIOMART format
hits[,4] <- paste(hits[,1], hits[,2], hits[,3], sep=":")
hitsfname <- paste(fname, ".biomart", sep="")
write.table(hits[,4], hitsfname, sep="\t", row.names=F, col.names=F, quote=F)











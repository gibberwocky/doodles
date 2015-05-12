setwd("/home/dave/Copy/HoneyBee/Analyses/selection/ehh")

fname <- "D-JFM-OTH-NC_007075.3.xpehh.out"
dat <- read.table(fname, header=T, stringsAsFactors=F)

# Standardize
dat[,8] <- (dat[,7]-mean(dat[,7]))/sd(dat[,7])

# Add chromosome and re-label vields
dat[,9] <- 6
names(dat)[9] <- "chr"
names(dat)[1] <- "start"


# Pick a value to plot
names(dat)[8] <- "Val"
plot_x_label <- "Z(XP-EHH)"
chromosome <- 0 # Set to 0 plots chromosomes 1-16
save2pdf <- 0

# Quantile
pc <- 0.1 # 0.1%
zHi_1 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)
pc <- 1 # 1%
zHi_2 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)

# ==============================================================================
# Generate plot
# ==============================================================================

# Sort data and get max and min Z_Hp values
tmp <- dat
if(chromosome > 0) tmp <- dat[which(dat$chr==chromosome),]
tmp <- subset(tmp[order(tmp$chr, tmp$start), ], )
ymax <- (max(tmp$Val, na.rm=T)) + (0.5*max(tmp$Val, na.rm=T))
ymin <- (min(tmp$Val, na.rm=T)) - (0.5*max(tmp$Val, na.rm=T))

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
if(save2pdf==1){
  pdfname <- paste(fname, "-", chromosome, ".pdf", sep="")
  pdf(pdfname, height=4, width=9.5)
}

# White background
bkgr="#FFFFFF"
frgr="#000000"
pcolz1=rgb(0.0,1.0,0.0)
pcolz2=rgb(0.0,0.8,0.0)
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
    cex=0.7, cex.sub=1, cex.main=1, cex.axis=1, cex.lab=1.5, lwd=1.5)
plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
     xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = xlabel, 
     ylab = plot_x_label, las = 1, pch = 20)
axis(1, at=ticks, labels=labs, lty=0)

# Points
col <- rep(rgb(0.2,0.2,0.2), nchr)
col[seq(1,nchr,2)] <- rgb(0.5,0.5,0.5)
#tmp[which(tmp$chr %% 2 == 0),3] <- rgb(0.2,0.2,0.2)
icol <- 1
for (i in unique(tmp$index)) {
  with(tmp[tmp$index == unique(tmp$index)[i], ], points(pos, Val, col = col[icol], pch = 20))
  icol = icol + 1
}
with(tmp[which(tmp$Val > zHi_2),], points(pos, Val, col = pcolz2, pch = 20))
with(tmp[which(tmp$Val > zHi_1),], points(pos, Val, col = pcolz1, pch = 20))
with(tmp[which(tmp$Val < 0- zHi_2),], points(pos, Val, col = pcolz2, pch = 20))
with(tmp[which(tmp$Val < 0- zHi_1),], points(pos, Val, col = pcolz1, pch = 20))

# Significance thresholds
abline(h=zHi_1, col=rgb(0.0,1.0,0.0))
abline(h=zHi_2, col=rgb(0.0,0.8,0.0))
abline(h=0-zHi_1, col=rgb(0.0,1.0,0.0))
abline(h=0-zHi_2, col=rgb(0.0,0.8,0.0))
title (main=fname)

# Known loci
ltype=1
lwidth=2
  # ROYAL JELLY
  # Dop2      LG15:6784429..6844358
  # hex7		  LG2:12690182..12694059
  # SsRbeta		LG10:13320..15501
  hex7 <- tmp[which(tmp$chr==2)[1],]$pos + 12690182
  SsRbeta <- tmp[which(tmp$chr==10)[1],]$pos + 13320
  Dop2 <- tmp[which(tmp$chr==15)[1],]$pos + 6784429
#  abline(v=hex7, col=rgb(0.5,0.5,1.0), lty=ltype, lwd=lwidth)
#  abline(v=SsRbeta, col=rgb(0.5,0.5,1.0), lty=ltype, lwd=lwidth)
#  abline(v=Dop2, col=rgb(0.5,0.5,1.0), lty=ltype, lwd=lwidth)

# MRJP genes (all between MRJPL and MRJPR)
MRJPL <-   tmp[which(tmp$chr==11)[1],]$pos + 2555000
abline(v=MRJPL, col=rgb(0.5,0.5,1.0), lty=ltype, lwd=lwidth)
MRJPR <-   tmp[which(tmp$chr==11)[1],]$pos + 2637000
abline(v=MRJPR, col=rgb(0.5,0.5,1.0), lty=ltype, lwd=lwidth)


  # FORAGING
  # For       LG13:8976929..9056381
  For <-  tmp[which(tmp$chr==13)[1],]$pos + 8976929
  abline(v=For, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
  # QTL from Page et al. 2012 - Complex pleiotropy...
  # stsD8 (pln1) QTL  LG13:3.4Mb
  stsD8 <-  tmp[which(tmp$chr==13)[1],]$pos + 3400000
  abline(v=stsD8, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
  # E15 (pln1) QTL  LG13:7.1Mb
  E15 <-  tmp[which(tmp$chr==13)[1],]$pos + 7100000
  abline(v=E15, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
  # W5 (pln2) QTL  LG1:16.9Mb
  W5 <-  tmp[which(tmp$chr==1)[1],]$pos + 16900000
  abline(v=W5, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
  # F8 (pln2) QTL  LG1:19.7Mb
  F8 <-  tmp[which(tmp$chr==1)[1],]$pos + 19700000
  abline(v=F8, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)
  # stsQ4 (pln3) QTL  LG1:9.2Mb
  stsQ4 <-  tmp[which(tmp$chr==1)[1],]$pos + 9200000
  abline(v=stsQ4, col=rgb(1.0,0.5,0.0), lty=ltype, lwd=lwidth)

box()

# Write to pdf
if(save2pdf==1) dev.off()


#pdfunite Fst-10k-25step.out.pdf Fst-20k-25step.out.pdf Fst-40k-25step.out.pdf \
#Fst-60k-25step.out.pdf Fst-80k-25step.out.pdf Fst-100k-25step.out.pdf \
#Fst-120k-25step.out.pdf Fst-140k-25step.out.pdf Fst-160k-25step.out.pdf \
#Fst.pdf





# ==============================================================================
# Write candidate regions to file
# ==============================================================================

#Subset to hits
hits <- tmp[which(tmp$Val>zHi_2),1:3]
  
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


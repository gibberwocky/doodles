library(stats)
library(pracma)

# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/selection")

# Integrated Haplotype Homozygosity standardized score
fname.ihs.jfm <- "iHS/OTH-16.ihs.out.100bins.norm" 
dat.ihs.jfm <- read.table(fname.ihs.jfm, header=F, stringsAsFactors=F)
# <locus> <pos> <'1' freq? <ihh1> <ihh0> <unstandardised iHS> <standardsised iHS> <is |iHS| >2 >
names(dat.ihs.jfm) <- c("SNP", "locus", "start", "frq1", "ihh1", "ihh0", "unstd_iHS", "Val", "iHS_above_2")


# Plot prep
bkgr="#FFFFFF"
frgr="#000000"
pcolz1=rgb(0.0,0.0,1.0)
pcolz2=rgb(0.8,0.8,1.0)
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
    cex=0.3, cex.sub=0.7, cex.main=0.7, cex.axis=0.7, cex.lab=1, lwd=1)
vrow <- 1
vcol <- 1
par(mfrow=c(vrow,vcol), oma = c(5.5,2.5,1,1), mar = c(1.5,0.5,1,0))
plot_x_label <- ""

# data prep
dat <- dat.ihs.jfm
tmp <- dat
tmp <- subset(tmp[order(tmp$start), ], )
ymax <- (max(tmp$Val, na.rm=T)) + (0.2*max(tmp$Val, na.rm=T))
ymin <- (min(tmp$Val, na.rm=T)) - (0.2*max(tmp$Val, na.rm=T))
xmin <- 0
xmax <- max(dat.ihs.jfm$start)

# Quantile
pc <- 0.01 # 0.01%
zHi_1 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)
pc <- 0.1 # 0.1%
zHi_2 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)

# draw plot
plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
     xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = NA, 
     ylab = plot_x_label, las = 1, pch = 19)

nchr <- 16
col <- rep(rgb(0.6,0.6,0.6), nchr)
col[seq(1,nchr,2)] <- rgb(0.8,0.8,0.8)
icol <- 1    

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

title(main="Integrated Haplotype Score, iHS (JFM)", line = 0.4)

title(xlab = "Position", outer=T)
axis(1, col=frgr, at=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))],
   labels=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))]) 



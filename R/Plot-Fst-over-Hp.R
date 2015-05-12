# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/selection")

# Pooled heterozygosity data JFM
fname.hp.jfm <- "JFM/1-Hp-100k-25step.out"    # A = 38.8456
dat.hp.jfm <- read.table(fname.hp.jfm, header=T, stringsAsFactors=F)
dat.hp.jfm <- dat.hp.jfm[which(dat.hp.jfm$n_SNPs>9),]
# Re-calculate Z_Hp across all autosomes rather than per chromosome as in file
dat.hp.jfm[which(dat.hp.jfm$chr < 17),11] <- 
  (dat.hp.jfm$Hp - mean(dat.hp.jfm[which(dat.hp.jfm$chr < 17),]$Hp, na.rm=T) ) / sd(dat.hp.jfm[which(dat.hp.jfm$chr < 17),]$Hp, na.rm=T)
names(dat.hp.jfm)[11] <- "Val"

# Pooled heterozygosity data OTH
fname.hp.oth <- "OTH/2-Hp-100k-25step.out"    # A = 38.8456
dat.hp.oth <- read.table(fname.hp.oth, header=T, stringsAsFactors=F)
dat.hp.oth <- dat.hp.oth[which(dat.hp.oth$n_SNPs>9),]
# Re-calculate Z_Hp across all autosomes rather than per chromosome as in file
dat.hp.oth[which(dat.hp.oth$chr < 17),11] <- 
  (dat.hp.oth$Hp - mean(dat.hp.oth[which(dat.hp.oth$chr < 17),]$Hp, na.rm=T) ) / 
  sd(dat.hp.oth[which(dat.hp.oth$chr < 17),]$Hp, na.rm=T)
names(dat.hp.oth)[11] <- "Val"

# Fst data
fname.fst <- "Fst-10k-25step.out"
dat.fst <- read.table(fname.fst, header=T, stringsAsFactors=F)
dat.fst <- dat.fst[which(dat.fst$n_SNPs>9),]
dat.fst[which(dat.fst$chr < 17),9] <- dat.fst$Fst_mean
names(dat.fst)[9] <- "Val"
plot_x_label <- "Fst"

# VCFT-Fst data
fname.vcftfst <- "VCFT-fst/S2-T4-JFM-OTH.windowed.weir.fst"
dat.vcftfst <- read.table(fname.vcftfst, header=T, stringsAsFactors=F)
dat.vcftfst <- dat.vcftfst[which(dat.vcftfst$N_VARIANTS>9),]
names(dat.vcftfst)[1] <- "chr"
names(dat.vcftfst)[2] <- "start"
names(dat.vcftfst)[5] <- "Val"

# XP-EHH data
fname.xpehh <- "/home/dave/Copy/HoneyBee/Analyses/selection/ehh/D-JFM-OTH-NC_007075.3.xpehh.out"
dat.xpehh <- read.table(fname.xpehh, header=T, stringsAsFactors=F)
dat.xpehh[,8] <- (dat.xpehh[,7]-mean(dat.xpehh[,7]))/sd(dat.xpehh[,7])
names(dat.xpehh)[1] <- "start"
names(dat.xpehh)[8] <- "Val"


# Add chromosome and re-label vields


n <- 6

#for(n in seq(1,16))
#{
  
chromosome <- n # Set to 0 plots chromosomes 1-16
savetopdf <- 0

# Update dat.xpehh chr field
dat.xpehh[,9] <- n
names(dat.xpehh)[9] <- "chr"


# Plot parameters
if(savetopdf==1){
  pdfname <- paste(chromosome, ".pdf", sep="")
  pdf(pdfname, height=9, width=6)
  png(pdfname)
}
bkgr="#FFFFFF"
frgr="#000000"
pcolz1=rgb(0.0,1.0,0.0)
pcolz2=rgb(0.0,0.8,0.0)
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
    cex=0.5, cex.sub=0.8, cex.main=0.7, cex.axis=0.8, cex.lab=1, lwd=1)
vrow <- 5
vcol <- 1
par(mfrow=c(vrow,vcol), oma = c(5.5,2.5,1,1), mar = c(1,0.5,0,0))



# ======================================================================================================
# PLOT
# ======================================================================================================

for(i in seq(1,5)){
  
if(i==1) dat <- dat.fst
if(i==2) dat <- dat.vcftfst
if(i==3) dat <- dat.hp.jfm
if(i==4) dat <- dat.hp.oth
if(i==5) dat <- dat.xpehh

# Quantile
pc <- 0.1 # 0.1%
zHi_1 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)
pc <- 1 # 1%
zHi_2 <- round(quantile(dat$Val, p=1-(pc/100), na.rm=T),2)

# Sort data and get max and min Z_Hp values
tmp <- dat[which(dat$chr==chromosome),]
tmp <- subset(tmp[order(tmp$chr, tmp$start), ], )
ymax <- (max(tmp$Val, na.rm=T)) + (0.2*max(tmp$Val, na.rm=T))
ymin <- (min(tmp$Val, na.rm=T)) - (0.2*max(tmp$Val, na.rm=T))
xmin <- 0
xmax <- max(dat.fst[which(dat.fst$chr==chromosome),]$start, 
            dat.hp.jfm[which(dat.hp.jfm$chr==chromosome),]$start, 
            dat.hp.oth[which(dat.hp.oth$chr==chromosome),]$start)

# Plot
plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
     xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = NA, 
     ylab = plot_x_label, las = 1, pch = 19)

# Points
nchr <- 17
col <- rep(rgb(0.2,0.2,0.2), nchr)
col[seq(1,nchr,2)] <- rgb(0.5,0.5,0.5)
icol <- 1
with(tmp, points(x=start, y=Val, col = col[icol], pch = 20))
with(tmp[which(tmp$Val > zHi_2),], points(x=start, y=Val, col = pcolz2, pch = 20))
with(tmp[which(tmp$Val > zHi_1),], points(x=start, y=Val, col = pcolz1, pch = 20))
with(tmp[which(tmp$Val < 0-zHi_2),], points(x=start, y=Val, col = pcolz2, pch = 20))
with(tmp[which(tmp$Val < 0-zHi_1),], points(x=start, y=Val, col = pcolz1, pch = 20))

# Significance thresholds
abline(h=zHi_1, col=rgb(0.0,1.0,0.0))
abline(h=zHi_2, col=rgb(0.0,0.8,0.0))
abline(h=0-zHi_1, col=rgb(0.0,1.0,0.0))
abline(h=0-zHi_2, col=rgb(0.0,0.8,0.0))


if(i==1) title (main="Fst (R)", line = -0.6)
if(i==2) title (main="Fst (VCFtools)", line = -0.6)
if(i==3) title (main="Pooled heterozygosity (JFM)", line = -0.6)
if(i==4) title (main="Pooled heterozygosity (OTH)", line = -0.6)
if(i==5) {
  title(main="XP-EHH (JFM | OTH)", line = -0.6)
  title(xlab = "Position", outer=T)
  axis(1, col=frgr, at=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))],
       labels=tmp$start[seq(1, length(tmp$start), by=(length(tmp$start)/20))]) 
}

# Box-it
box()

}

if(savetopdf==1) dev.off()

#}


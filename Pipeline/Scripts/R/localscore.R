#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(IRanges)
library(stats)
library(pracma)
library(fdrtool)
library(mnormt)
#fsttmp <- ( dat.fst$Val - mean(dat.fst$Val, na.rm=T) ) / sd(dat.fst$Val, na.rm=T)
#test <- fdrtool(fsttmp)

# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/selection")


# Pooled heterozygosity data OTH
fname.hp.oth <- "Hp/5kb/JFMOTH_preImp_chrs16_qc_beagle-OTH.vcf.hp"
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
dat.li.oth[,12] <- 2*pnorm(-abs(dat.li.oth$ZHp)) # Converts Z score to P value which already has uniform distribution
names(dat.li.oth)[12] <- "UZHp"
#pc <- 0.01 # Set significance threshold for Linldey process at 0.01%
#OTHlizHi_1 <- round(quantile(dat.li.oth$UZHp, p=1-(pc/100), na.rm=T),2)
OTHlizHi_1 <- quantile(dat.li.oth$UZHp, .95)
dat.li.oth[,13] <- -log10(dat.li.oth$UZHp) - OTHlizHi_1 # Create holder field for storing results of Lindley process
# Perform Linldey process on each chromosome independently to avoid overlapping calculations

# Create vector to store Rho values
rho.chr <- NULL

# Forward direction
for(i in unique(dat.li.oth$chr)){
  scores <- -log10(dat.li.oth[which(dat.li.oth$chr==i),]$UZHp) - OTHlizHi_1  
  h=rep(NA,length(scores))
  h[1]=0   
  for (n in 2:length(scores))  h[n]=max(0,h[n-1]+scores[n])
  dat.li.oth[which(dat.li.oth$chr==i),13] <- h
  zz <- dat.li.oth[which(dat.li.oth$chr==i),]$UZHp
  rho.chr <- c(rho.chr, cor(dat.li.oth[which(dat.li.oth$chr==i),]$UZHp[-1], dat.li.oth[which(dat.li.oth$chr==i),]$UZHp[-length(scores)]))
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
#dat.li.oth[,15] <- rowMeans(dat.li.oth[,13:14])
#names(dat.li.oth)[15] <- "Val"

# Standardise local score values (may or may not be sensible thing to do)
#dat.li.oth[,13] <- (dat.li.oth[,13] - mean(dat.li.oth[,13], na.rm=T) ) /  sd(dat.li.oth[,13], na.rm=T)
#dat.li.oth[,14] <- (dat.li.oth[,14] - mean(dat.li.oth[,14], na.rm=T) ) /  sd(dat.li.oth[,14], na.rm=T)
#dat.li.oth[,15] <- (dat.li.oth[,15] - mean(dat.li.oth[,15], na.rm=T) ) /  sd(dat.li.oth[,15], na.rm=T)



### simulation of correlated uniform variables from correlated Gaussian ###
fscore = function(x){return(-log10(x)-OTHlizHi_1)}
runifgausscor = function ( n, rho)
{
  out = outX= X= rnorm(1)
  for (i in 1:n) {
    X=rnorm(1, mean=rho*X, sd=sqrt(1-rho*rho))
    outX = c(outX, X)}
  #out = ecdf(out)(out)
  out=pnorm(outX)
  return(out)
}

pval = runifgausscor( 1000, rho=mean(rho.chr) )
hist(pval) #ok
plot(acf(pval)) #autocorrelation
cor(pval[1:999],pval[2:1000]) # ok
plot(pval) 



################  H0  #####################
nsim=10000  #number of runs
rho=mean(rho.chr)     #autocorrelation
nloc=5000   #nimber of SNP (length of sequence)

max.score=scoreloc.cont=rho.obs=rep(NA,nsim)
for (isim in 1:nsim)
{
  pval=runifgausscor(nloc, rho=rho)
  rho.obs[isim]=cor(pval[-nloc],pval[-1])
  score.cont=fscore(pval)
  h=rep(NA,nloc)
  h[1]=max(0,score.cont[1])
  #h[1]=0   #ce n'est pas exactement ca, mais sinon ca coince plus loin ...
  for (i in 2:nloc)  h[i]=max(0,h[i-1]+score.cont[i])
  scoreloc.cont[isim]=max(h)
}

rho.theor=pmnorm( c(0,0), c(0,0), matrix(c(1,rho/2,rho/2,1),ncol=2) )
rho.theor=12*rho.theor-3
rho.theor

hist(rho.obs)
abline(v=rho,col="green",lwd=2)
abline(v=rho.theor,col="blue",lwd=2)
abline(v=mean(rho.obs),col="black",lwd=2,lty=2)
quantile(scoreloc.cont,0.95) 








# Compute observed autocorrelation
scores <- dat.li.oth[,13]

a <- log(length(scores)) - 7.6 - 6.86*rho
b <- -0.49 + 1.51 * rho
alpha <- 0.05
thres = ( log(-log(1-alpha)) - a ) / b

plot(dat.li.oth[,13])
abline(h=thres, col="red", lwd=2)







OTH <- dat.hp.oth[which(dat.hp.oth$chr==6),]

png("Notes/Hist-Hp.png", height=6, width=6, units="in", res=300)
par(mar=c(4,4,3,3), oma=c(2,2,2,2), cex=0.8, cex.main=1)
xlim = c(quantile(OTH$Hp,0), quantile(OTH$Hp,1))
h <- hist(OTH$Hp, breaks=25, freq=F, xlim=xlim, col="grey", xlab=expression(bold(paste(italic(H[P])))), main="")
title (main = expression(bold(paste("Pooled heterozygosity", sep=""))), line = 0.4)
dev.off()

png("Notes/Hist-ZHp.png", height=6, width=6, units="in", res=300)
par(mar=c(4,4,3,3), oma=c(2,2,2,2), cex=0.8, cex.main=1)
xlim = c(quantile(OTH$Val,0), quantile(OTH$Val,1))
h <- hist(OTH$Val, breaks=25, freq=F, xlim=xlim, col="grey", xlab=expression(bold(paste(italic(Z(H[P]))))), main="")
title (main = expression(bold(paste("Standardised pooled heterozygosity", sep=""))), line = 0.4)
dev.off()


png("Notes/Hp.png", height=6, width=9, units="in", res=300)
col <- rep(rgb(0.85,0.85,0.85), 16)
col[seq(1,16,2)] <- rgb(0.8,0.8,0.8)
icol <- 1    
zHi_1 <- 3.719
xmin <- 1
xmax <- max(OTH$start)
ymin <- min(0-zHi_1-1, OTH$Val-1)
ymax <- max(zHi_1+1, OTH$Val+1)
plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
     xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = NA, 
     ylab = "", las = 1, pch = 19)
with(OTH, points(x=start, y=Val,  col = col[icol], pch=20, cex=0.5))

with(OTH[which(OTH$Val > zHi_1),], points(x=start, y=Val, col = "blue", pch = 20))
with(OTH[which(OTH$Val < 0-zHi_1),], points(x=start, y=Val, col = "blue", pch = 20))
abline(h=zHi_1, col="blue")
abline(h=0-zHi_1, col="blue")
abline(h=round(mean(OTH$Val, na.rm=T),3), col=rgb(0.3,0.3,0.3), lty=2, lwd=1)    

title(xlab = paste("Chromosome 6, position", sep=""), outer=F)
axis(1, col="black", at=OTH$start[seq(1, length(OTH$start), by=(length(OTH$start)/20))],
     labels=OTH$start[seq(1, length(OTH$start), by=(length(OTH$start)/20))]) 
title(ylab = expression(italic(Z(H[P]))), outer=F, line= 2)
title (main= expression(bold(paste(italic(Z(H[P])), sep=""))), line = 0.4)
dev.off()




OTH <- dat.li.oth[which(dat.li.oth$chr==6),]

png("Notes/Hist-UZHp.png", height=6, width=6, units="in", res=300)
par(mar=c(4,4,3,3), oma=c(2,2,2,2), cex=0.8, cex.main=1)
xlim = c(quantile(OTH$UZHp,0), quantile(OTH$UZHp,1))
h <- hist(OTH$UZHp, breaks=25, freq=F, xlim=xlim, col="grey", xlab=expression(bold(paste(2 * ~ pnorm(-abs(italic(Z(H[P]))))))), main="")
title (main = expression(bold(paste("Reverse standardised pooled heterozygosity", sep=""))), line = 0.4)
dev.off()

png("Notes/LocalScore.png", height=6, width=9, units="in", res=300)
col <- rep(rgb(0.85,0.85,0.85), 16)
col[seq(1,16,2)] <- rgb(0.8,0.8,0.8)
icol <- 1    
zHi_1 <- 3.719
xmin <- 1
xmax <- max(OTH$start)
ymin <- min(OTH$Val)
ymax <- max(c(OTH$fVal, OTH$rVal, OTH$Val))
plot(NULL, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
     xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = NA, 
     ylab = "", las = 1, pch = 19)
with(OTH, points(x=start, y=OTH$rVal, type="l", lwd=2, lty=3, col=rgb(1,0.5,0.5,0.5)))
with(OTH, points(x=start, y=OTH$fVal, type="l", lwd=2, lty=3, col=rgb(0.5,1,0.5,0.5)))  
with(OTH, points(x=start, y=Val,  col = col[icol], pch=20, cex=0.5))
if( is.null(nrow(OTH[which(OTH$Val< zHi_1),]$Val))==F) OTH[which(OTH$Val< zHi_1),]$Val <- NA
with(OTH[which(OTH$Val >= zHi_1),], points(x=start, y=Val, pch=20, col = "blue", cex=0.5))
abline(h=zHi_1, col="blue", lwd=1)        
title(xlab = paste("Chromosome 6, position", sep=""), outer=F)
axis(1, col="black", at=OTH$start[seq(1, length(OTH$start), by=(length(OTH$start)/20))],
     labels=OTH$start[seq(1, length(OTH$start), by=(length(OTH$start)/20))]) 
title(ylab = expression(italic(Z(H[L]))), outer=F, line= 2)
title (main= expression(bold(paste(italic(Z(H[P])), " local score", sep=""))), line = 0.4)
dev.off()






setwd("/home/dave/Copy/HoneyBee/Analyses/JFM-OTH-Apis/Admixture-LD03-BP5K")
tmp <- read.table("JFM-OTH-BR-Apis.fam", header=F, stringsAsFactors=F)[,1]
ids <- NULL
for(i in 1:length(tmp)){
  ids <- c(ids,unlist(strsplit(tmp, split="-")[i])[1])
}

K2 <- read.table("JFM-OTH-BR-Apis.2.Q", header=F, stringsAsFactors=F)
K3 <- read.table("JFM-OTH-BR-Apis.3.Q", header=F, stringsAsFactors=F)
K4 <- read.table("JFM-OTH-BR-Apis.4.Q", header=F, stringsAsFactors=F)
K5 <- read.table("JFM-OTH-BR-Apis.5.Q", header=F, stringsAsFactors=F)

# Re-order the SRA populations, makes it easier to interpret ancestry
tmp <- cbind(K2, K3, K4, K5)
tmp2 <- tmp[3:55,]
rownames(tmp2) <- ids[3:55]
tmp <- tmp[-(3:55),]
rownames(tmp) <- c(ids[1:2], ids[56:length(ids)])

# Sort ancestral populations into mitotype groupings
tmp.1 <- tmp2[which(tmp2[,14]>0.5),]
tmp.1 <- tmp.1[order(tmp.1[,14], decreasing=T),]

tmp.2 <- tmp2[which(tmp2[,13]>0.5),]
tmp.2 <- tmp.2[order(tmp.2[,13], decreasing=T),]

tmp.3 <- tmp2[which(tmp2[,12]>0.5),]
tmp.3 <- tmp.3[order(tmp.3[,12], decreasing=T),]

tmp.4 <- tmp2[which(tmp2[,11]>0.5),]
tmp.4 <- tmp.4[order(tmp.4[,11], decreasing=T),]

# Stick all of data back together
tmp2 <- rbind(tmp.1, tmp.2, tmp.3, tmp.4, tmp)

K2 <- tmp2[,1:2]
K3 <- tmp2[,3:5]
K4 <- tmp2[,6:9]
K5 <- tmp2[,10:14]
ids <- rownames(tmp2)
  
# For dark background
#frgr <- "#FFFFFF" # white
#bkgr <- "#333333" # dkgrey
#pop1 <- rgb(1.0,0.5,0.0) # orange
#pop2 <- rgb(0.5,0.5,1.0) # blue
#pop3 <- rgb(1.0,1.0,0.0) # yellow
#pop4 <- rgb(0.0,0.7,0.3) # green
#pop5 <- rgb(0.0,1.0,1.0) # cyan

# For white background
frgr <- "#000000" # black
bkgr <- "#FFFFFF" # white
pop1 <- "#ff0000" # unknown
pop2 <- "#02b7cf" # O  Turqoise
pop3 <- "#fdc701" # C  Orange
pop4 <- "#df71ab" # A  Pink
pop5 <- "#000000" # M  Black

# Hex colour codes
# A #df71ab
# M #000000
# M #ff0000 (iberiensis)
# C #fdc701
# O #bebebe (syriaca)
# O #02b7cf (anatolica)


# K5
pdf("admixture.pdf", height=6, width=12)
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
    cex=0.5, cex.sub=1, cex.main=1, cex.axis=1, cex.lab=1.5, lwd=1.5)
vrow <- 2
vcol <- 1
par(mfrow=c(vrow,vcol),
    oma = c(5.5,2.5,1,1) ,
    mar = c(1,0.5,0,0))
par(cex=1, cex.main=0.8, cex.lab=0.8, cex.axis=0.6, cex.sub=0.7, las=3)
#barplot(t(as.matrix(K2)), col=c(pop1,pop2),
#        xlab=NA, ylab="Ancestry", border=1)
#barplot(t(as.matrix(K3)), col=c(pop3,pop2,pop1),
#        xlab=NA, ylab="Ancestry", border=1)
barplot(t(as.matrix(K4)), col=c(pop2,pop1,pop3,pop4),
        xlab=NA, ylab="Ancestry", border=1)
barplot(t(as.matrix(K5)), col=c(pop1,pop2,pop5,pop4,pop3),
        xlab=NA, ylab="Ancestry", border=1, names=ids)
title (ylab="Q-value", outer=T, line=1.5)

# K1, n=9
X <- c(0.5, 10.5)
axis(1, at=X, col="red", line=4.2,
     tick=T, labels=rep("",2), lwd=2, lwd.ticks=0) 
mtext("C", 1, at=(X[2]-X[1])/2, las=1, line=4.3)

# K2, n=25
X <- c(11, 40.5)
axis(1, at=X, col="red", line=4.2, 
     tick=T, labels=rep("",2), lwd=2, lwd.ticks=0) 
mtext("A", 1, at=(X[2]-X[1])/2 + X[1], las=1, line=4.3)

# K3, n=9
X <- c(41, 51.5)
axis(1, at=X, col="red", line=4.2,
     tick=T, labels=rep("",2), lwd=2, lwd.ticks=0) 
mtext("M", 1, at=(X[2]-X[1])/2 + X[1], las=1, line=4.3)

# K4, n=10
X <- c(52, 63.5)
axis(1, at=X, col="red", line=4.2,
     tick=T, labels=rep("",2), lwd=2, lwd.ticks=0) 
mtext("O", 1, at=(X[2]-X[1])/2 + X[1], las=1, line=4.3)

dev.off()

#



# K3 no labels
pdf("admixture.pdf", height=6, width=12)
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
    cex=0.5, cex.sub=1, cex.main=1, cex.axis=1, cex.lab=1.5, lwd=1.5)
vrow <- 3
vcol <- 1
par(mfrow=c(vrow,vcol),
    oma = c(0.2,0.2,0.2,0.2) ,
    mar = c(0.5,0.5,0,0))
par(cex=1, cex.main=0.8, cex.lab=0.8, cex.axis=0.6, cex.sub=0.7, las=3)
barplot(t(as.matrix(K2)), col=c(pop1,pop2), ,
        xlab=NA, ylab=NA, axes=F, border=1)
barplot(t(as.matrix(K3)), col=c(pop3,pop2,pop1),
        xlab=NA, ylab=NA, axes=F, border=1)
barplot(t(as.matrix(K4)), col=c(pop2,pop1,pop4,pop3),
        xlab=NA, ylab=NA, axes=F, border=1)
dev.off()

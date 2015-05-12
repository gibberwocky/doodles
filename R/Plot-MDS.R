setwd("/home/dave/Copy/HoneyBee/Analyses/JFM-OTH-Apis/Admixture-LD03-BP5K")
mds <- read.table("plink.mds", header=T, stringsAsFactors=F)
tmp <- read.table("JFM-OTH-BR-Apis.fam", header=F, stringsAsFactors=F)[,1]
ids <- NULL
for(i in 1:length(tmp)){
  ids <- c(ids,unlist(strsplit(tmp, split="-")[i])[1])
}
rownames(mds) <- ids


K2 <- read.table("JFM-OTH-BR-Apis.2.Q", header=F, stringsAsFactors=F)
K3 <- read.table("JFM-OTH-BR-Apis.3.Q", header=F, stringsAsFactors=F)
K4 <- read.table("JFM-OTH-BR-Apis.4.Q", header=F, stringsAsFactors=F)
K5 <- read.table("JFM-OTH-BR-Apis.5.Q", header=F, stringsAsFactors=F)

# For white background
frgr <- "#000000" # black
bkgr <- "#FFFFFF" # white
pop1 <- rgb(1.0,0.0,0.5) # purple
pop2 <- rgb(1.0,0.4,0.0) # dark orange
pop3 <- rgb(1.0,0.7,0.1) # light orange
pop4 <- rgb(1.0,1.0,0.1) # orange-yellow
pop5 <- rgb(0.0,0.7,0.2) # green

# Re-order the SRA populations, makes it easier to interpret ancestry
tmp <- cbind(K2, K3, K4, K5)
tmp2 <- tmp[3:55,]
rownames(tmp2) <- ids[3:55]
tmp <- tmp[-(3:55),]
rownames(tmp) <- c(ids[1:2], ids[56:length(ids)])

# Sort ancestral populations into mitotype groupings
tmp.1 <- tmp2[which(tmp2[,14]>0.5),]
tmp.2 <- tmp2[which(tmp2[,13]>0.5),]
tmp.3 <- tmp2[which(tmp2[,12]>0.5),]
tmp.4 <- tmp2[which(tmp2[,11]>0.5),]

# Colour blobs
mds[rownames(tmp.1),]$SOL <- pop3
mds[rownames(tmp.2),]$SOL <- pop4
mds[rownames(tmp.3),]$SOL <- pop5
mds[rownames(tmp.4),]$SOL <- pop2
mds[1:2,]$SOL <- "#000000"
mds[56:85,]$SOL <- pop3
mds[86:nrow(mds),]$SOL <- pop1

# MDS PLOT
pdf("mds.pdf", height=12, width=6)

par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, 
    cex=0.8, cex.sub=0.8, cex.main=0.8, cex.axis=0.8, cex.lab=0.8, lwd=1.5)
vrow <- 3
vcol <- 1
par(mfrow=c(vrow,vcol),
    oma = c(1, 1, 1, 1) ,
    mar = c(3,3,1,1))
par(cex=0.8, cex.main=0.8, cex.lab=0.8, cex.axis=0.6, cex.sub=0.7, las=3)

gap <- 2
plot(x=mds$C1, y=mds$C2, col=mds$SOL, pch=20, xlab="n", ylab="n")
text(mds$C1, mds$C2, rownames(mds), pos=3, cex=0.3)
title(ylab="CV2", outer=F, line=gap, cex=0.8)
title(xlab="CV1", outer=F, line=gap, cex=0.8)

plot(x=mds$C1, y=mds$C3, col=mds$SOL, pch=20, xlab="n", ylab="n")
text(mds$C1, mds$C3, rownames(mds), pos=3, cex=0.3)
title(ylab="CV3", outer=F, line=gap, cex=0.8)
title(xlab="CV1", outer=F, line=gap, cex=0.8)

plot(x=mds$C2, y=mds$C3, col=mds$SOL, pch=20, xlab="n", ylab="n")
text(mds$C2, mds$C3, rownames(mds), pos=3, cex=0.3)
title(ylab="CV3", outer=F, line=gap, cex=0.8)
title(xlab="CV2", outer=F, line=gap, cex=0.8)

dev.off()




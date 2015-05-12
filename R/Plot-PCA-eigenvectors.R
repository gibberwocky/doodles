library(ape)
library(cluster) 
library(plotrix)
library(car)

# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/admixture")

# Import eigen-vectors
eigs <- read.table("Admixture-dip-prune.eigenvec", header=F, strings=F, sep=" ")
rownames(eigs) <- eigs[,1]
eigs <- eigs[,3:ncol(eigs)]
colnames(eigs) <- paste("PCA", seq(1,ncol(eigs)), sep="-")
km <- kmeans(eigs[,1:2], centers=4, iter.max=1000)
C_rgb <- "#fdc701"
A_rgb <- "#df71ab"
M_rgb <- "#000000"
O_rgb <- "#02b7cf"
pdf("pca.pdf", height=6, width=6)
px <- 1
py <- 2
par(cex=0.8)
plot(eigs[,c(px,py)], type="n")
  points(eigs[names(km$cluster[which(km$cluster==1)]),c(px,py)], col=A_rgb, pch=20)
  points(eigs[names(km$cluster[which(km$cluster==2 & substr(names(km$cluster),1,1)=="S")]),c(px,py)], col=C_rgb, pch=20)
  points(eigs[names(km$cluster[which(km$cluster==2 & substr(names(km$cluster),1,1)!="S")]),c(px,py)], col=C_rgb, pch=17)
  points(eigs[names(km$cluster[which(km$cluster==3)]),c(px,py)], col=O_rgb, pch=20)
  points(eigs[names(km$cluster[which(km$cluster==4)]),c(px,py)], col=M_rgb, pch=20)

# For PCA1,PCA2
  cluster_id <- c("A", "C", "O", "M")
  text(km$centers, labels=cluster_id, pos=4)
dev.off()


# Cladogram
pdf("cladogram.pdf", height=6, width=6)
bkgr="#FFFFFF"
frgr="#000000"
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, cex=0.6, cex.sub=0.7, cex.main=0.7, cex.axis=0.7, cex.lab=0.7, lwd=1)
phy <- as.phylo(hclust(dist(eigs[,1:2])))
plot(phy,type="cladogram", label.offset=0.002, font=1,
     edge.color=c(rep(A_rgb,22), rep(M_rgb,17), A_rgb, rep(O_rgb,19), rep(C_rgb,77)) )
dev.off()








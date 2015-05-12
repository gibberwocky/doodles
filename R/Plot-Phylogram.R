library(ape)
library(cluster) 
library("GenABEL")

# Set working directory
setwd("/home/dave/Copy/HoneyBee/Analyses/admixture/R")

# Import PED data
ped <- "Admixture-dip-prune-R.ped"
map <- "Admixture-dip-prune-R.map"
pheno <- "Admixture-dip-prune-R.pheno"

convert.snp.ped(pedfile=ped, mapfile=map, out="gen.raw", format = "premakeped", traits = 1, strand = "u", bcast = 10000000, wslash=F, mapHasHeaderLine=F)
raw <- load.gwaa.data(phe=pheno, gen = "gen.raw", force = T, sort=T)


raw@phdata[which(substr(raw@phdata$id,1,1)=="S"),]
#phy <- as.phylo(hclust(dist(as.numeric(raw[which(substr(raw@phdata$id,1,1)=="S"),]@gtdata))))
phy <- as.phylo(hclust(dist(as.numeric(raw@gtdata))))

# Colours for branches
C_rgb <- "#fdc701"
A_rgb <- "#df71ab"
M_rgb <- "#000000"
O_rgb <- "#02b7cf"
# n Branches = number of edge.length
# All data
#branches <- c(rep("#000000",length(phy$edge.length)))
#branches[1:78] <- C_rgb
#branches[79:95] <- M_rgb
#branches[96:115] <- O_rgb
#branches[116:136] <- A_rgb

# [which(substr(raw@phdata$id,1,1)=="S"),]
branches <- c(rep("#000000",length(phy$edge.length)))
#branches[1:18] <- C_rgb
#branches[19:35] <- M_rgb
#branches[36:55] <- O_rgb
#branches[56:76] <- A_rgb
branches[60:136] <- C_rgb
#branches[1:18] <- M_rgb
branches[19:38] <- O_rgb
branches[c(1,39:59)] <- A_rgb

# Plot
pdf("Harpur-unrooted.pdf", height=6, width=6)
bkgr="#FFFFFF"
frgr="#000000"
par(bg=bkgr, fg=frgr, col.lab=frgr, col.axis=frgr, cex=0.6, cex.sub=0.7, cex.main=0.7, cex.axis=0.7, cex.lab=0.7, lwd=1)
plot(phy,type="unrooted", label.offset=1.5, font=1, edge.color=branches, lab4ut="axial" )
dev.off()








# PCA
pca <- prcomp(as.numeric(raw@gtdata), weight="freq")
km <- kmeans(pca$x[,1:2], centers=4, iter.max=1000, nstart=1000)
px <- 1
py <- 2

pdf("pca.pdf", height=6, width=6)
par(cex=0.8)
plot(pca$x[,c(px,py)], type="n")
points(pca$x[names(km$cluster[which(km$cluster==2)]),c(px,py)], col=O_rgb, pch=20)
points(pca$x[names(km$cluster[which(km$cluster==1 & substr(names(km$cluster),1,1)=="S")]),c(px,py)], col=C_rgb, pch=20)
points(pca$x[names(km$cluster[which(km$cluster==1 & substr(names(km$cluster),1,1)!="S")]),c(px,py)], col=C_rgb, pch=2)
points(pca$x[names(km$cluster[which(km$cluster==3)]),c(px,py)], col=M_rgb, pch=20)
points(pca$x[names(km$cluster[which(km$cluster==4)]),c(px,py)], col=A_rgb, pch=20)
# For PCA1,PCA2
cluster_id <- c("C", "O", "M", "A")
text(km$centers, labels=cluster_id, pos=4)
dev.off()





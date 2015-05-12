#setwd("/home/dave/Copy/HoneyBee/Analyses/JFM-OTH-Apis/Admixture-LD03-BP5K")
#setwd("/home/dave/Copy/HoneyBee/Analyses/JFM-OTH-Apis/test-LD03-BP5K")
#setwd("/home/dave/Copy/HoneyBee/Analyses/JFM-OTH-Apis/test-geno0-LD03-BP5K")
#setwd("/home/dave/Copy/HoneyBee/Analyses/JFM-OTH-Apis/test-hom")
setwd("/home/dave/Copy/HoneyBee/Analyses/JFM-OTH-Apis/AMEL-HARP-MISC")
#tmp <- read.table("JFM-OTH-BR-Apis.fam", header=F, stringsAsFactors=F)[,1]
tmp <- read.table("Amel-test_A.fam", header=F, stringsAsFactors=F)[,1]
ids <- NULL
for(i in 1:length(tmp)){
  ids <- c(ids,unlist(strsplit(tmp, split="-")[i])[1])
}

#K2 <- read.table("JFM-OTH-BR-Apis.2.Q", header=F, stringsAsFactors=F)
#K3 <- read.table("JFM-OTH-BR-Apis.3.Q", header=F, stringsAsFactors=F)
#K4 <- read.table("JFM-OTH-BR-Apis.4.Q", header=F, stringsAsFactors=F)
#K5 <- read.table("JFM-OTH-BR-Apis.5.Q", header=F, stringsAsFactors=F)

K2 <- read.table("Amel-test_A.2.Q", header=F, stringsAsFactors=F)
K3 <- read.table("Amel-test_A.3.Q", header=F, stringsAsFactors=F)
K4 <- read.table("Amel-test_A.4.Q", header=F, stringsAsFactors=F)
K5 <- read.table("Amel-test_A.5.Q", header=F, stringsAsFactors=F)
K6 <- read.table("Amel-test_A.6.Q", header=F, stringsAsFactors=F)

# K5 Fst table from Admixture
#     OTH   O     M     A	
# O   0.426	
# M   0.367	0.437	
# A   0.284	0.293	0.284	
# C   0.118 0.381 0.340 0.255	


# Re-order the SRA populations, makes it easier to interpret ancestry
tmp <- cbind(K2, K3, K4, K5, K6)
tmp2 <- tmp[c(1,2:52),]
rownames(tmp2) <- ids[c(1,2:52)]
#tmp <- tmp[-(c(1,2:52)),]
#rownames(tmp) <- c(ids[2:3], ids[57:length(ids)])

# Sort ancestral bees
tmp2[,(ncol(tmp2)+1)] <- 0
for(x in 1:nrow(tmp2)) tmp2[x,ncol(tmp2)] <- which(tmp2[x,]==max(tmp2[x,6:9]))
tmp2 <- tmp2[order(tmp2[,21],tmp2[,6],tmp2[,7],tmp2[,8],tmp2[,9],decreasing=T),]
colnames(tmp2) <- paste("V", seq(1,21), sep="")
#colnames(tmp) <- paste("V", seq(1,20), sep="")
pops <- tmp2[,21]
#tmp2 <- rbind(tmp2[,1:20], tmp)

K2 <- tmp2[,1:2]
K3 <- tmp2[,3:5]
K4 <- tmp2[,6:9]
K5 <- tmp2[,10:14]
K6 <- tmp2[,15:20]
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
pop3 <- "#ff0000" # unknown
pop2 <- "#02b7cf" # O  Turqoise
pop1 <- "#fdc701" # C  Orange
pop4 <- "#df71ab" # A  Pink
pop5 <- "#000000" # M  Black
pop6 <- "blue"   # Unknown

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
vrow <- 3
vcol <- 1
par(mfrow=c(vrow,vcol),
    oma = c(5.5,2.5,1,1) ,
    mar = c(1,0.5,0,0))
par(cex=1, cex.main=0.8, cex.lab=0.5, cex.axis=0.5, cex.sub=0.5, las=3)
#barplot(t(as.matrix(K2)), col=c(pop1,pop2),
#        xlab=NA, ylab="Ancestry", border=1)
#barplot(t(as.matrix(K3)), col=c(pop3,pop2,pop1),
#        xlab=NA, ylab="Ancestry", border=1)
barplot(t(as.matrix(K4)), col=c(pop4,pop5,pop1, pop2),
        xlab=NA, ylab="Ancestry", border=NA, space=0)


# Annotations
inds <- c(rev(as.vector(table(pops))))
labs <- c("O", "C", "M", "A")

for(i in 1:length(inds)){
  if(i==1) X <- c(0.2, inds[1]-0.2)
  if(i>1) X <- c(sum(inds[1:(i-1)])+0.2, sum(inds[1:i])-0.2)
  axis(1, at=X, col="red", line=0.2,
       tick=T, labels=rep("",2), lwd=2, lwd.ticks=0)   
  rect(X[1]-0.2, 0, X[2]+0.2, 1)
  mtext(labs[i], 1, at=X[1] + (X[2]-X[1])/2, las=1, line=0.05, cex=0.7)
}

barplot(t(as.matrix(K5)), col=c(pop2,pop1,pop5,pop4,pop3),
        xlab=NA, ylab="Ancestry", border=NA, space=0)

barplot(t(as.matrix(K6)), col=c(pop5,pop6,pop1,pop4,pop2,pop3),
        xlab=NA, ylab="Ancestry", border=NA, space=0, names=ids)

title (ylab="Q-value", outer=T, line=1.5)




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

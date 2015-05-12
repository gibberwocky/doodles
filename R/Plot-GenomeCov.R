fname <- "/home/dave/Ubuntu One/HoneyBee/Notes/JFM11/JFM11_GGCTAC_L001_realn.bam.dcov"

# Flagstat
tmp <- read.delim(fname, sep="\t", stringsAsFactors=F, header=F)
data <- tmp[which(tmp[,1]=="genome"),]
data[nrow(data),2] <- paste("≥ ", data[nrow(data),2], sep="")

# Plot
tiff(paste(fname, ".tif", sep=""), compression="lzw")
options("scipen"=100, "digits"=4)
par(cex=1, cex.main=0.9, cex.lab=0.9, cex.axis=0.9, cex.sub=0.9)

ymax <- round(max(as.numeric(data[,5]))*1.1,1)
ymin <- 0

barplot(as.numeric(data[,5]), ylim=c(ymin, ymax), names.arg=data[,2], col=rainbow(nrow(data)))
grid(nx=NA, ny=10, lty=1, col=rgb(0.4,0.4,0.8), lwd=0.25)
barplot(as.numeric(data[,5]), ylim=c(ymin, ymax), names.arg=data[,2], col=rainbow(nrow(data)), add=T)

title (xlab="Depth of coverage", line=3)
title (ylab="Percentage of genome", line=3)
title (main="JFM11_GGCTAC_L001", line=3, cex.main = 0.9)
box()

dev.off()



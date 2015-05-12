#!/usr/bin/env Rscript

args <- commandArgs(T)

# (1) sample name
# (2) path/sample_realn.flagstat

# Read and parse flagstat data
tmp <- read.delim(args[2], sep=" ", stringsAsFactors=F, header=F)
data <- data.frame(n=numeric(nrow(tmp)-1), t=character(nrow(tmp)-1), stringsAsFactors=F)
for(i in 1:nrow(data)){
	data[i,1] <- tmp[i,1]
	a <- paste(tmp[i,4:ncol(tmp)])
	if(i==nrow(data)) a <- c(a, tmp[nrow(tmp),1])
	data[i,2] <- paste(a, collapse=" ")
}
data <- data[c(1, 2, 3, 4, 7, 8, 5, 6, 9, 10, 11),]

# Plot
pdf(paste(args[2], ".pdf", sep=""), height=5, width=5)
options("scipen"=100, "digits"=4)
par(cex=0.9, cex.main=0.5, cex.lab=0.5, cex.axis=0.5, cex.sub=0.5)
buffer <- 1.1
ymax <- round(max(as.numeric(data[,1]))*1.1,1)
ymin <- 0
barplot(as.numeric(data[,1]), ylim=c(ymin, ymax), col=rainbow(nrow(data)))
legend("topright", data[,2], col = rainbow(nrow(data)), text.col = rgb(0.0,0.0,0.0), pch = rep(19, nrow(data)), pt.lwd=2, cex=0.5)
title(xlab="Metric", line=1)
title(ylab="Number of reads", line=2)
title(main=args[1])
box()
dev.off()




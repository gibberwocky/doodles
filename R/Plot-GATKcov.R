#!/usr/bin/env Rscript

args <- commandArgs(T)

# (1) sample name
# (2) path/sample_realn_GATKcov.sample_summary


# GATKcoverage summary
tmp <- read.table(args[2], sep="\t", stringsAsFactors=F, header=T, nrows=1, strip.white=T)
data <- t(tmp[,3:6])

# Plot
pdf(paste(args[2], ".pdf", sep=""), height=5, width=5)
options("scipen"=100, "digits"=4)
par(cex=0.9, cex.main=0.5, cex.lab=0.5, cex.axis=0.5, cex.sub=0.5)
buffer <- 1.1
ymax <- round(max(as.numeric(data[,1]))*1.1,1)
ymin <- 0
barplot(as.numeric(data[,1]), ylim=c(ymin, ymax), col=rainbow(nrow(data)))
legend("topright", rownames(data), col = rainbow(nrow(data)), text.col = rgb(0.0,0.0,0.0), pch = rep(19, nrow(data)), pt.lwd=2, cex=0.5)
title(main=paste("GATK DepthOfCoverage (", args[1], ")", sep=""))
title(xlab="Metric", line=1)
title(ylab="Depth of coverage", line=2)
box()
dev.off()



# Pct_Bases covered by X 
data <- t(tmp[,7:9])

# Plot
pdf(paste(args[2], "_base.pdf", sep=""), height=5, width=5)
options("scipen"=100, "digits"=4)
par(cex=0.9, cex.main=0.5, cex.lab=0.5, cex.axis=0.5, cex.sub=0.5)
buffer <- 1.1
ymax <- round(max(as.numeric(data[,1]))*1.1,1)
ymin <- 0
barplot(as.numeric(data[,1]), ylim=c(ymin, ymax), col=rainbow(nrow(data)))
legend("topright", rownames(data), col = rainbow(nrow(data)), text.col = rgb(0.0,0.0,0.0), pch = rep(19, nrow(data)), pt.lwd=2, cex=0.5)
title(main=paste("GATK, % bases covered by X number of reads (", args[1], ")", sep=""))
title(xlab="Dpeth of coverage", line=1)
title(ylab="Percentage of bases", line=2)
box()
dev.off()



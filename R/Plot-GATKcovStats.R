#!/usr/bin/env Rscript

args <- commandArgs(T)

# (1) sample name
# (2) path/sample_realn_GATKcov.sample_statistics


# GATKcoverage summary
tmp <- read.table(args[2], sep="\t", stringsAsFactors=F, header=T, strip.white=T)
data <- t(tmp[,2:50])
data <- rbind(data, sum(tmp[,51:ncol(tmp)]))
rownames(data) <- seq(1, 50)
ymax <- round(max(as.numeric(data[,1]))*1.1,1)
ymin <- 0
pdf(paste(args[2], "_interval.pdf", sep=""), height=5, width=5)
par(cex=0.9, cex.main=0.5, cex.lab=0.5, cex.axis=0.5, cex.sub=0.5)
barplot(t(data), ylim=c(ymin, ymax), col=rgb(0,0.5,1))
title(main=paste("GATK DepthOfCoverage (", args[1], ")", sep=""))
title(xlab="Depth of coverage", line=2)
title(ylab="Number of bases", line=2)
box()
dev.off()




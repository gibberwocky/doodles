#!/usr/bin/env Rscript

args <- commandArgs(T)

# (1) output path
# (2) sample name
# (3) path/sample_realn.bam.dcov
# (4) path/sample_gatk.sample_summmary

library(Hmisc)

# Create a copy of last 21 lines of file(3)
tmp <- paste(args[3], ".tmp", sep="")
var_sys	<- paste("tail --l 21 ", args[3], " > ", tmp, sep="")
system(var_sys)

# Read in copy of file(3), and file(4)
coverage_X <- read.table(tmp)[,c(2,5)]
coverage_mean <- read.table(args[4], header=F, skip=3)[3]

# Calculate poisson distribution and ymax for plot
pk <- dpois(seq(0,nrow(coverage_X)), as.numeric(coverage_mean))
ymax <- round(max(coverage_X[,2], pk, na.rm=T),2)+0.1

# Generate distribution plot of genome coverage
fname <- paste(args[1], "/", args[2], "/metrics/", args[2], "_Coverage_Poisson.pdf", sep="")
pdf(fname, height=5, width=5)
options("scipen"=100, "digits"=4)
par(cex=0.6)
mp<-barplot(coverage_X[,2], names.arg = coverage_X[,1], ylim=c(0, ymax), xlab="X depth of coverage", ylab="fraction of bases (Amel 4.5)", col=rgb(0,1,0,0.5))
lines(pk, lwd=2, lty=3, col=rgb(0.5,0,0,0.5))                  
title(main=args[2])
legend("topright", c(args[2], 
  as.expression(substitute(paste(lambda," = ",var),list(var=as.numeric(coverage_mean))) )),
  col=c(rgb(0,1,0,0.5), rgb(0.5,0,0,0.5)),
  text.col=rgb(0,0,0), pch=c(19,18), pt.lwd=2, cex=0.75)
box()
dev.off()


# Record P-value of correlation between Poisson and Empirical
fname <- paste(args[1], "/", args[2], "/metrics/", args[2], "_Coverage_Poisson.pval", sep="")
Pval <- rcorr(cbind(pk, coverage_X[,2]), type="pearson")
write(Pval$P[1,2], fname, append=F, sep="\t")




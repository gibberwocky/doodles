#!/usr/bin/env Rscript

args <- commandArgs(T)

# VennDiagram library
library("VennDiagram")

setwd("/home/dave/Copy/HoneyBee/Manuscripts/SeqApiPop/GeneSets")

alphaJFM <- as.vector(read.table("alphaJFM", header=F, strings=F)[,1])
alphaOTH <- as.vector(read.table("alphaOTH", header=F, strings=F)[,1])
FST <- as.vector(read.table("FST", header=F, strings=F)[,1])
HPJFMLoss <- as.vector(read.table("HPJFMLoss", header=F, strings=F)[,1])
HPJFMExec <- as.vector(read.table("HPJFMExec", header=F, strings=F)[,1])
HPOTHLoss <- as.vector(read.table("HPOTHLoss", header=F, strings=F)[,1])
XPEHHanc <- as.vector(read.table("XPEHHanc", header=F, strings=F)[,1])
XPEHHder <- as.vector(read.table("XPEHHder", header=F, strings=F)[,1])

# alpha, FST
venn.diagram(list(alphaJFM=alphaJFM, alphaOTH=alphaOTH, FST=FST),
	fill=c("red", "blue", "green"),alpha=c(0.5, 0.5, 0.5), cex=0.5, cat.cex=0.6,
	cat.fontfamily="sans", cat.fontface=2, fontfamily="sans", fontface=1,
	lty=1, euler.d=TRUE, scaled=TRUE, reverse=TRUE, compression='lzw', margin=0.1,
	resolution=300, height=6, width=6, units='in', filename="alpha.tiff");

# Fst - XPEHH
venn.diagram(list(FST=FST, XPEHHanc=XPEHHanc, XPEHHder=XPEHHder),
             fill=c("red", "blue", "green"),alpha=c(0.5, 0.5, 0.5), cex=0.5, cat.cex=0.6,
             cat.fontfamily="sans", cat.fontface=2, fontfamily="sans", fontface=1,
             lty=1, euler.d=TRUE, scaled=TRUE, reverse=TRUE, compression='lzw', margin=0.1,
             resolution=300, height=6, width=6, units='in', filename="fst-xpehh.tiff");

# HP
venn.diagram(list(HPJFMLoss=HPJFMLoss, HPJFMExec=HPJFMExec, HPOTHLoss=HPOTHLoss),
             fill=c("red", "blue", "green"),alpha=c(0.5, 0.5, 0.5), cex=0.5, cat.cex=0.6,
             cat.fontfamily="sans", cat.fontface=2, fontfamily="sans", fontface=1,
             lty=1, euler.d=TRUE, scaled=TRUE, reverse=TRUE, compression='lzw', margin=0.1,
             resolution=300, height=6, width=6, units='in', filename="hp.tiff");




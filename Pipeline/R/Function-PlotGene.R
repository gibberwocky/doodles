


plotGene <- function(target, gene, hap1.N, hap1.S, bp1, bp2, pop1.name="JFM", pop2.name="OTH", pop1, pop2,
                     fname="", type = c("png", "pdf"), height=4, width=4, units="in", res=600, gap=7,
                     ref=rgb(0.95,0.95,0.95), alt.N=rgb(0.4,0.4,1.0), alt.S=rgb(0.4,0.9,0.4), het=rgb(1.0,0.4,0.4), 
                     gene.col.exon=rgb(0.7,0.7,1.0), gene.col.cds=rgb(0.2,0.2,1.0), gene.col.start=rgb(0.4,1.0,0.4),
                     gene.col.stop=rgb(1.0,0.4,0.4), bg="#FFFFFF", fg="#000000") {

  ref=rgb(0.95,0.95,0.95); alt.N=rgb(0.4,0.4,1.0); alt.S=rgb(0.4,0.9,0.4); het=rgb(1.0,0.4,0.4)
  gene.col.exon=rgb(0.7,0.7,1.0); gene.col.cds=rgb(0.2,0.2,1.0); gene.col.start=rgb(0.4,1.0,0.4)
  gene.col.stop=rgb(1.0,0.4,0.4); bg="#FFFFFF"; fg="#000000"; gap=7

  # Plot details
  if(fname != ""){
    if (type=="png") png(fname, height=height, width=width, units=units, res=res)
    if (type=="pdf") pdf(fname, height=height, width=width)
  }
  
  
  # Plot y size
  ysize <- ncol(hap1.N)+gap
  par(bg=bg, fg=fg, col.lab=fg, col.axis=fg, cex=0.3, cex.sub=0.4, cex.main=0.4, cex.axis=0.3, cex.lab=0.4, lwd=1)
  par(mfrow=c(1,1), oma = c(3,3,1,1), mar = c(0.5,0.5,0.5,0.5))
  plot(NULL, yaxt = "n", xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", xlim = c(bp1, bp2), 
       ylim = c(0, ysize), xlab = NA, ylab = NA, las = 1, pch = 19)
  N.allele.color <- c(ref, het, alt.N)
  S.allele.color <- c(ref, het, alt.S)
  
  # Gene
  rect(bp1, ysize-0.9, bp2, ysize-1.2, col=fg, border=NA)
  mtext(target, 2, at=ysize+0.5, las=1, line=0.5, cex=0.5, adj=1)
  
  # Exons
  tmp <- lapply(rownames(gene[which(gene$feature=="exon"),]), function(i) rect(gene[i,]$start, ysize-0.5, gene[i,]$end, ysize-1.5, col=gene.col.exon, border=fg))
  # CDS
  tmp <- lapply(rownames(gene[which(gene$feature=="CDS"),]), function(i) rect(gene[i,]$start, ysize-0.5, gene[i,]$end, ysize-1.5, col=gene.col.cds, border=fg))
  # Start codon
  tmp <- lapply(rownames(gene[which(gene$feature=="start_codon"),]), function(i) rect(gene[i,]$start, ysize-0.5, gene[i,]$end, ysize-1.5, col=gene.col.start, border=fg))
  # Stop codon
  tmp <- lapply(rownames(gene[which(gene$feature=="stop_codon"),]), function(i) rect(gene[i,]$start, ysize-0.5, gene[i,]$end, ysize-1.5, col=gene.col.stop, border=fg))
  
  
  # Replace chromosome with flag for non-synonymous / synonymous for reference once merged
  if(nrow(hap1.N)>0) hap1.N$chr <- 1
  if(nrow(hap1.S)>0) hap1.S$chr <- 0
  
  # SNP legs
  hap <- (rbind(hap1.N, hap1.S))
  hap <- t(hap[order(hap$pos),])
  snps.spaced <- seq(bp1+10, bp2-10, (bp2-bp1)/(nrow(hap1.N)+nrow(hap1.S)+2) )
  snps.spaced <- snps.spaced[2:(length(snps.spaced)-1)]
  
  for(y in 1:ncol(hap)){
    segments(hap[2,y], ysize-2, hap[2,y], ysize-3, col=fg, lty=3)
    segments(hap[2,y], ysize-3, snps.spaced[y], ysize-5, col=fg, lty=3)
    segments(snps.spaced[y], ysize-5, snps.spaced[y], ysize-7, col=fg, lty=3)
  }
  
  # Alleles Pop 1
  pop.hap <- data.frame(hap[pop1,])
  pop1.size <- nrow(pop.hap)
  if(ncol(pop.hap)>1){
    pop.dist <- hclust(dist(pop.hap))
    pop.hap <- pop.hap[c(pop.dist$order),]
    colnames(pop.hap) <- colnames(hap)
    pop.hap <- rbind(hap[1:2,], pop.hap)
  } else {
    pop.hap[1:(nrow(pop.hap)+2),1] <- c(hap[1,1], hap[2,1], pop.hap[1:nrow(pop.hap),1])
    rownames(pop.hap) <- c("chr", "pos", rownames(pop.hap)[1:(nrow(pop.hap)-2)])
  }
  for (y in 1:ncol(pop.hap)){
    for (x in 3:nrow(pop.hap)) { 
      if (pop.hap[1,y]==0) { points(x=snps.spaced[y], y=(ysize-6)-x, col=S.allele.color[pop.hap[x,y]+1], pch=15, cex=0.7)}
      if (pop.hap[1,y]==1) { points(x=snps.spaced[y], y=(ysize-6)-x, col=N.allele.color[pop.hap[x,y]+1], pch=15, cex=0.7)}
    }
  }
  mtext(pop1.name, 2, at=(ysize-6)-(nrow(pop.hap)/2), las=1, line=0.5, cex=0.5, adj=1)
  
  # Divider
  segments(bp1, (ysize-6)-pop1.size-3, bp2, (ysize-6)-pop1.size-3, lty=3, col=fg)
  
  # Alleles Pop 1
  pop.hap <- data.frame(hap[pop2,])
  if(ncol(pop.hap)>1){
    pop.dist <- hclust(dist(pop.hap))
    pop.hap <- pop.hap[c(pop.dist$order),]
    colnames(pop.hap) <- colnames(hap)
    pop.hap <- rbind(hap[1:2,], pop.hap)
  } else {
    pop.hap[1:(nrow(pop.hap)+2),1] <- c(hap[1,1], hap[2,1], pop.hap[1:nrow(pop.hap),1])
    rownames(pop.hap) <- c("chr", "pos", rownames(pop.hap)[1:(nrow(pop.hap)-2)])
  }
  for (y in 1:ncol(pop.hap)){
    for (x in 3:nrow(pop.hap)) { 
      if (pop.hap[1,y]==0) { points(x=snps.spaced[y], y=(ysize-7)-x-pop1.size, col=S.allele.color[pop.hap[x,y]+1], pch=15, cex=0.7)}
      if (pop.hap[1,y]==1) { points(x=snps.spaced[y], y=(ysize-7)-x-pop1.size, col=N.allele.color[pop.hap[x,y]+1], pch=15, cex=0.7)}
    }
  }
  mtext(pop2.name, 2, at=(ysize-6-pop1.size)-(nrow(pop.hap)/2), las=1, line=0.5, cex=0.5, adj=1)
  
  if(fname != "") dev.off()
  
}
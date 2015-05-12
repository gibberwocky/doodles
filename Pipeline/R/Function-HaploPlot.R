


haploplot <- function(target, gene, hap1, hap2, hap.N, hap.S, bp1, bp2, extend=0, pop1.name="JFM", pop2.name="OTH", pop1, pop2,
                     fname="", type = c("png", "pdf"), height=4, width=4, units="in", res=600, gap=7,
                     ref=rgb(0.95,0.95,0.95), alt=rgb(0.7,0.7,0.7), het=rgb(1.0,0.4,0.4), 
                     alt.N=rgb(0.2,0.2,0.8), alt.S=rgb(0.4,0.9,0.4), gene.col.exon=rgb(0.7,0.7,1.0),
                     gene.col.cds=rgb(0.2,0.2,1.0), gene.col.start=rgb(0.4,1.0,0.4), gene.col.stop=rgb(1.0,0.4,0.4), 
                     bg="#FFFFFF", fg="#000000") {
  
  
  
  #target=dat.hits[hit,1]; hap1=jfm.haps; hap2=oth.haps; pop1.name="JFM"; pop2.name="OTH"; pop1=JFM; pop2=OTH
  #type="png"; height=175; width=175; units="mm"; res=600; gap=7; extend=500
  #ref=rgb(0.95,0.95,0.95); alt=rgb(0.7,0.7,0.7); het=rgb(1.0,0.4,0.4)
  #alt.N=rgb(0.2,0.2,0.8); alt.S=rgb(0.4,0.9,0.4); gene.col.exon=rgb(0.7,0.7,1.0)
  #gene.col.cds=rgb(0.2,0.2,1.0); gene.col.start=rgb(0.4,1.0,0.4); gene.col.stop=rgb(1.0,0.4,0.4)
  #bg="#FFFFFF";fg="#000000"
  
  
  
  pos1 <- max(0,bp1-extend)
  pos2 <- min(bp2+extend, as.numeric(colnames(hap1)[ncol(hap1)]))
  hapTmp <- data.frame(hap1[,which(as.numeric(colnames(hap1))>=pos1 & as.numeric(colnames(hap1))<=pos2)],check.names=F)
  
  if(ncol(hapTmp)>0){
    
    # Plot details
    if(fname != ""){
      if (type=="png") png(fname, height=height, width=width, units=units, res=res)
      if (type=="pdf") pdf(fname, height=height, width=width)
    }
    
    # Plot y size
    ysize <- ncol(hap1.N)+gap
    par(bg=bg, fg=fg, col.lab=fg, col.axis=fg, cex=0.3, cex.sub=0.4, cex.main=0.4, cex.axis=0.3, cex.lab=0.4, lwd=1)
    par(mfrow=c(1,1), oma = c(3,3,1,1), mar = c(0.5,0.5,0.5,0.5))
    plot(NULL, yaxt = "n", xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", xlim = c(pos1, pos2), 
         ylim = c(0, ysize), xlab = NA, ylab = NA, las = 1, pch = 19)
    allele.color <- c(ref, het, alt)
    N.allele.color <- c(ref, het, alt.N)
    S.allele.color <- c(ref, het, alt.S)
    
    # Strand
    if(gene$strand[1]=="-") segments(pos1, ysize-2, pos2, ysize-2, col="darkgreen", lty=3, lwd=1)
    if(gene$strand[1]=="+") segments(pos1, ysize-0.1, pos2, ysize-0.1, col="red", lty=3, lwd=1)
    
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
    snps.spaced <- seq(pos1+10, pos2-10, (pos2-pos1)/(ncol(hapTmp)+2) )
    snps.spaced <- snps.spaced[2:(length(snps.spaced)-1)]
    if(length(snps.spaced)>1){
      for(y in 1:ncol(hapTmp)){
        segments(as.numeric(colnames(hapTmp)[y]), ysize-3, as.numeric(colnames(hapTmp)[y]), ysize-4, col=fg, lty=3)
        segments(as.numeric(colnames(hapTmp)[y]), ysize-4, snps.spaced[y], ysize-5, col=fg, lty=3)
        segments(snps.spaced[y], ysize-5, snps.spaced[y], ysize-6, col=fg, lty=3)
      }
    } else { 
      colnames(hapTmp)[1] <- snps.spaced[1]
      segments(snps.spaced[1],ysize-3, as.numeric(colnames(hapTmp)[1]), ysize-6, col=fg, lty=3) 
    }
    #rm(hapTmp)
    
    # JFM haplotypes
    dat <- hap1[,which(as.numeric(colnames(hap1))>=pos1 & as.numeric(colnames(hap1))<=pos2)]
    dat.t <- as.data.frame(t(dat))
    dat.dist <- hclust(dist(dat))
    dat.t <- dat.t[,c(dat.dist$order)]
    dat.t[,(ncol(dat.t)+1)] <- rownames(dat.t)
    names(dat.t)[ncol(dat.t)] <- "pos"
    for(i in 1:(ncol(dat.t)-1)) {
      for(x in 1:nrow(dat.t)) { 
        colAllele <- allele.color[dat.t[x,i]+1]
        if(as.numeric(dat.t[x,]$pos) %in% hap.N$pos) colAllele <- N.allele.color[dat.t[x,i]+1]
        if(as.numeric(dat.t[x,]$pos) %in% hap.S$pos) colAllele <- S.allele.color[dat.t[x,i]+1]
        points(snps.spaced[x], (ysize-gap)-i/2, col = colAllele, pch=15, cex=0.5) 
      }
    }
    mtext(pop1.name, 2, at=(ysize-gap)-(nrow(hap1)/4), las=1, line=0.5, cex=0.5, adj=1)
    
  
    # Horizontal divider
    segments(pos1, (ysize-gap)-(i/2)-0.75, pos2, (ysize-gap)-(i/2)-0.75, lwd=1, col=fg, lty=3)
    
    
    # OTH haplotypes
    dat <- hap2[,which(as.numeric(colnames(hap2))>=pos1 & as.numeric(colnames(hap2))<=pos2)]
    dat.t <- as.data.frame(t(dat))
    dat.dist <- hclust(dist(dat))
    dat.t <- dat.t[,c(dat.dist$order)]
    dat.t[,(ncol(dat.t)+1)] <- rownames(dat.t)
    names(dat.t)[ncol(dat.t)] <- "pos"
    for(i in 1:(ncol(dat.t)-1)) {
      for(x in 1:nrow(dat.t)) { 
        colAllele <- allele.color[dat.t[x,i]+1]
        if(as.numeric(dat.t[x,]$pos) %in% hap.N$pos) colAllele <- N.allele.color[dat.t[x,i]+1]
        if(as.numeric(dat.t[x,]$pos) %in% hap.S$pos) colAllele <- S.allele.color[dat.t[x,i]+1]
        points(snps.spaced[x], (ysize-gap)-(i/2)-(nrow(hap1)/2)-1, col = colAllele, pch=15, cex=0.5) 
      }
    }
    mtext(pop2.name, 2, at=(ysize-gap)-(nrow(hap2)/4)-(i/2), las=1, line=0.5, cex=0.5, adj=1)
    
    if(fname != "") dev.off()
    return(1)
  } else {
    return(0)
  }
  
}
  




plotHaps <- function(genes, hap1, hap2, start, end, extend=0, gap=6, pop1="JFM", pop2="OTH",
                     fname="", type = c("png", "pdf"), height=4, width=12, res=600,
                     ref=rgb(0.4,1.0,0.4), alt=rgb(0.4,0.4,1.0), het=rgb(1.0,0.4,0.4), bg="#FFFFFF", fg="#000000") {
  
  
  # Configure positions
  extend <- extend
  bp1 <- start - extend
  bp2 <- end + extend      
  map <- as.numeric(colnames(hap1))
  tmp <- which( map >= bp1 & map <= bp2)
  
  if(length(tmp)>0){
    pos1 <- as.numeric(min(tmp))
    pos2 <- as.numeric(max(tmp))  
      
    # Plot details
    if(fname != ""){
      if (type=="png") png(fname, height=height, width=width, units="in", res=res)
      if (type=="pdf") pdf(fname, height=height, width=width)
    }
    
    # draw plot
    allele.color <- c(ref, het, alt)
    par(bg=bg, fg=fg, col.lab=fg, col.axis=fg, cex=0.3, cex.sub=0.4, cex.main=0.4, cex.axis=0.3, cex.lab=0.4, lwd=1)
    par(mfrow=c(1,1), oma = c(2,2,1,1), mar = c(1.5,0.5,1,0))
  
    plot(NULL, yaxt = "n", xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", xlim = c(bp1, bp2), 
         ylim = c(0, nrow(hap1)+gap), xlab = NA, ylab = NA, las = 1, pch = 19)
    
    # JFM haplotypes
    dat <- hap1[,pos1:pos2]
    dat.t <- as.data.frame(t(dat))
    dat.dist <- hclust(dist(dat))
    dat.t <- dat.t[,c(dat.dist$order)]
    dat.t[,(ncol(dat.t)+1)] <- rownames(dat.t)
    names(dat.t)[ncol(dat.t)] <- "pos"
    for(i in 1:(ncol(dat.t)-1)) {
      for(x in 1:nrow(dat.t)) { points(dat.t[x,]$pos, i/2, col = allele.color[dat.t[x,i]+1],pch=15, cex=0.5) }
    }
    mtext(pop1, 2, at=8, las=1, line=0.5, cex=0.5)
      
    # OTH haplotypes
    dat <- hap2[,pos1:pos2]
    dat.t <- as.data.frame(t(dat))
    dat.dist <- hclust(dist(dat))
    dat.t <- dat.t[,c(dat.dist$order)]
    dat.t[,(ncol(dat.t)+1)] <- rownames(dat.t)
    names(dat.t)[ncol(dat.t)] <- "pos"
    for(i in 1:(ncol(dat.t)-1)) {
      for(x in 1:nrow(dat.t)) { points(dat.t[x,]$pos, 17+i/2, col = allele.color[dat.t[x,i]+1],pch=15, cex=0.5) }
    }
    mtext(pop2, 2, at=23, las=1, line=0.5, cex=0.5)
        
    # Annotation
    text <- paste(dat.type, ", Haplotypes: ", bp2-bp1+1, " SNPs", sep="")
    title(main=text, line = 0.4)
    xtext <- paste("Chr ", chromosome, ": position", sep="")
    title(xlab = xtext, outer=T, line=-0.6)
    axis(1, col=fg, at=seq(bp1, bp2, by=(bp2-bp1)), labels=F, line=0, tcl=-0.2) 
    axis(1, col=fg, at=seq(bp1, bp2, by=(bp2-bp1)), tick=F, line=-1, tcl=-0.2) 
    axis(1, col=fg, at=dat.t$pos, labels=names(dat), tick=F, line=-0.75, tcl=-0.2) 
    axis(1, col=fg, at=dat.t$pos, labels=NA, tick=T, line=0, tcl=-0.2) 
    
    
    
    # Gene horizontal line break
    for.colour <- rgb(0.0,0.8,0.0)
    rev.colour <- rgb(0.8,0.0,0.0)
    mtext("Genes", 1, at=bp1, line=1.5, cex=0.4, col="black")
    mtext("F", 1, at=bp1, line=2, cex=0.3, col=for.colour)
    mtext("R", 1, at=bp1, line=2.4, cex=0.3, col=rev.colour)
    axis(1, at=c(bp1, bp2), col="blue", line=2.5, tick=T, labels=rep("",2), lwd=1, lwd.ticks=0)   
    axis(1, at=c(bp1, bp2), col="blue", line=2.9, tick=T, labels=rep("",2), lty=2, lwd=0.5, lwd.ticks=0)       
    
    # Genes which overlap with the plotted interval
    int.genes <- genes[which(genes$chr == chromosome & genes$end >= bp1 & genes$start <= bp2),]
    # Forward genes
    for.genes <- int.genes[which(int.genes$strand=="+"),]
    for.exons <- int.genes[which(int.genes$strand=="+" & int.genes$feature=="CDS"),]
    for.codns <- int.genes[which(int.genes$strand=="+" & (int.genes$feature=="start_codon" | int.genes$feature=="stop_codon")),]
    if (nrow(for.genes)>0){
      for(i in 1:nrow(for.genes)) {
        for.min <- bp1
        for.max <- bp2
        for.att <- strsplit(for.genes[i,]$attribute, split="; ")
        for.att <- strsplit(for.att[[1]][1], split=" ")
        for.att <- for.att[[1]][length(for.att[[1]])]
        for.min <- min(for.genes[which(grepl(for.att, for.genes$attribute)!=F), 4:5])
        if (sum(for.genes[which(grepl(for.att, for.genes$attribute)!=F),]$feature %in% "start_codon")==0) for.min <- bp1
        for.max <- max(for.genes[which(grepl(for.att, for.genes$attribute)!=F), 4:5])
        if (sum(for.genes[which(grepl(for.att, for.genes$attribute)!=F),]$feature %in% "stop_codon")==0) for.max <- bp2
        rect(for.min, 16.7, for.max, 16.8, col=for.colour, border=NA)
        mtext(for.att, 1, at=(for.min + ((for.max-for.min)/2)), cex=0.3, line=2, col=for.colour)
      }
      if(nrow(for.exons)>0) for(i in 1:nrow(for.exons)) {
        rect(for.exons[i,]$start, 16.6, for.exons[i,]$end, 16.9, col=for.colour, border=NA)
      }
      if(nrow(for.codns)>0) for(i in 1:nrow(for.codns)) {
        rect(for.codns[i,]$start, 16.6, for.codns[i,]$end, 16.9, col="yellow", border=T, density=-1)
      }
    }
    
    # Reverse genes
    rev.genes <- int.genes[which(int.genes$strand=="-"),]
    rev.exons <- int.genes[which(int.genes$strand=="-" & int.genes$feature=="CDS"),]
    rev.codns <- int.genes[which(int.genes$strand=="-" & (int.genes$feature=="start_codon" | int.genes$feature=="stop_codon")),]
    if (nrow(rev.genes)>0){
      for(i in 1:nrow(rev.genes)) {
        rev.min <- bp1
        rev.max <- bp2
        rev.att <- strsplit(rev.genes[i,]$attribute, split="; ")
        rev.att <- strsplit(rev.att[[1]][1], split=" ")
        rev.att <- rev.att[[1]][length(rev.att[[1]])]
        rev.max <- min(rev.genes[which(grepl(rev.att, rev.genes$attribute)!=F), 4:5])
        if (sum(rev.genes[which(grepl(rev.att, rev.genes$attribute)!=F),]$feature %in% "start_codon")==0) rev.min <- bp2
        rev.min <- max(rev.genes[which(grepl(rev.att, rev.genes$attribute)!=F), 4:5])
        if (sum(rev.genes[which(grepl(rev.att, rev.genes$attribute)!=F),]$feature %in% "stop_codon")==0) rev.max <- bp1
        rect(rev.max, 16.3, rev.min, 16.4, col=rev.colour, border=NA)
        mtext(rev.att, 1, at=(rev.min + ((rev.max-rev.min)/2)), cex=0.3, line=2.4, col=rev.colour)
      }
      if(nrow(rev.exons)>0) for(i in 1:nrow(rev.exons)) {
        rect(rev.exons[i,]$start, 16.2, rev.exons[i,]$end, 16.5, col=rev.colour, border=NA)
      }
      if(nrow(rev.codns)>0) for(i in 1:nrow(rev.codns)) {
        rect(rev.codns[i,]$start, 16.2, rev.codns[i,]$end, 16.5, col="yellow", border=T, density=-1)
      }  
    }
    
    if(fname != "") dev.off()
    return(1)
  } else {
    return(0)
  }
}

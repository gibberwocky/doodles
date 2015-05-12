# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/selection/VEP")

#fname <- "fst.bed.clean"
#fname <- "hp-OTH.bed.clean"
#fname <- "hp-JFM.bed.clean"
fname <- "XPEHH.bed.clean"

# Read table of hits (this is after having merged the BED files)
dat.hits <- read.table(fname, sep="\t", header=F, strings=F)
names(dat.hits) <- c("chr", "start", "end", "n-sum", "method", "min", "max")

# Read table of GO terms (this is after VEP prediction, etc.)
dat.GOs <- read.table("mart_export-GOs.txt", sep="\t", header=T, strings=F)
names(dat.GOs) <- c("Gene", "chr", "start", "end", "GO", "Name", "Domain")

# ReiGO table
# Have to extend the search position by 5 Kb to capture ontologies which have been picked up from
# transcripts in addition to gene coordinates
dist <- 5000
revigo <- NULL
tmp.err <- NULL
x<-0
for(i in unique(dat.hits$chr)){
  cat("Chr:", i,"\n")
  tmp1 <- dat.hits[which(dat.hits$chr==i),]
  tmp2 <- dat.GOs[which(dat.GOs$chr==i),]
  tmp2 <- tmp2[-which(tmp2$GO==""),]
  
  for(n in 1:nrow(tmp2)){
    
    #tmp3 <- tmp2[which( (tmp2$start >= tmp1[n,]$start & tmp2$start <= tmp1[n,]$end) | 
    #                      (tmp2$end >= tmp1[n,]$start & tmp2$end <= tmp1[n,]$end)),]

    tmp3 <- tmp1[which( (tmp1$start >= (tmp2[n,]$start-dist) & tmp1$start <= (tmp2[n,]$end+dist)) | 
                          (tmp1$end >= (tmp2[n,]$start-dist) & tmp1$end <= (tmp2[n,]$end+dist)) |
                          (tmp1$start <= tmp2[n,]$start & tmp1$end >= tmp2[n,]$end)),]
    
    
    if(nrow(tmp3)>0) {
      x<-x+1
#      tmp3[,ncol(tmp3)+1] <- tmp1[n,]$min
#      tmp3[,ncol(tmp3)+1] <- tmp1[n,]$max
  tmp3[,ncol(tmp3)+1] <- tmp2[n,]$GO
tmp3[,ncol(tmp3)+1] <- tmp2[n,]$Gene
      revigo <- rbind(revigo, tmp3)
    }
    if(nrow(tmp3)==0){
      tmp.err <- rbind(tmp.err, tmp2[n,])
    }
  }
  
}
revigo$min <- round(revigo$min, 3)
revigo$max <- round(revigo$max, 3)
names(revigo)[ncol(revigo)-1] <- "GO"
names(revigo)[ncol(revigo)] <- "Gene"

# Record largest absolute distance from 0 for the value
revigo[,10] <- 0
names(revigo)[10] <- "Val"
for(i in 1:nrow(revigo)){
  revigo[i,]$Val <- max(abs(0-revigo[i,]$min), revigo[i,]$max)
}

# Extract unique GO terms and the max Val associated with them
revigo.agg <- merge(aggregate(Val ~ GO, revigo, max), revigo)

# Write revigo, and the first to columns of revigo.agg (for REViGO analysis)
write.table(revigo.agg[,1:2], paste(fname, ".rev-agg.txt", sep=""), sep="\t", append=F, quote=F, row.names=F, col.names=F)
write.table(revigo, paste(fname, ".rev", sep=""), sep="\t", append=F, quote=F, row.names=F, col.names=F)






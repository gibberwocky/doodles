# Set path
setwd("/home/dave/Copy/HoneyBee/Analyses/selection/DAF")
setwd("/work/dwragg/Analysis/selection/daf")

fname.daf <- "S2_JFM-OTH.daf"
dat <- read.table(fname.daf, header=T, stringsAsFactors=F)
names(dat) <- c("chr", "start", "AF1", "AF2", "DAF", "zDAF")

winL <- 5000
winO <- 1000

DAF <- NULL
for(i in unique(dat$chr))
{
  cat("CHR: ",i,"\n")
  # Calculates intervals based on parameters
  tmp <- dat[which(dat$chr==i),]
  chrL <- max(tmp$start)
  start <- as.integer(seq(1,chrL, by=winO))
  chr <- rep(i, length(start))
  stop <- as.integer(start + winL)
  uDAF <- rep(NA, length(start))
  bed <- NULL
  bed <- data.frame(chr, start, stop, uDAF, stringsAsFactors=F)
  bed[nrow(bed),ncol(bed)] <- chrL
  a <- lapply(as.integer(rownames(bed)), 
               function(n) mean(tmp[which(tmp$start >= bed[n,]$start & tmp$start <= bed[n,]$stop),]$DAF, na.rm=T)) 
  bed$uDAF <- unlist(a)  
  
#  for(n in 1:nrow(bed)) { 
#    bed[n,]$uDAF <- mean(tmp[which(tmp$start >= bed[n,]$start & tmp$start <= bed[n,]$stop),]$DAF, na.rm=T)  
#  }

  DAF <- rbind(DAF, bed)
}

DAF[,5] <- (DAF[,4]-mean(DAF[,4],na.rm=T)/sd(DAF[,4],na.rm=T))
names(DAF)[5] <- "zDAF"
fname <- paste(fname.daf, ".win", sep="")
write.table(DAF, fname, sep="\t", col.names=T, row.names=F,quote=F, append=F)



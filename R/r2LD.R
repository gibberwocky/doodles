# r2LD.R

# Read in .ld file output by Plink r2 function
setwd("~/work/Analysis/vcfs/ld")
pop <- 2
# Get list of chromosome files
chrs <- list.files()

# Configure bins
bins <- seq(0,1000000,1000)
#bins <- c(0, 50, 100, 200, 400, 800, 1200, 1600, 2000, 
#	2500, 3000, 3500, 4000, 4500, 5000)

# Dataframe
data <- as.data.frame(t(matrix(nrow=length(chrs),ncol=length(bins))))

# Iterate through each chromosome
for(i in 1:length(chrs))
{
	names(data)[i] <- chrs[i]
	LD <- read.table(chrs[i], header=T, stringsAsFactors=F)
	for(n in 1:length(bins))
	{
		if(n<length(bins)){
		    data[n,i] <- round(mean(LD[which( 
			((LD$BP_B - LD$BP_A) >= bins[n]) & 
			((LD$BP_B - LD$BP_A) < bins[n+1]) ),]$R2, na.rm=T),2)
		} else {
		    data[n,i] <- round(mean(LD[which( 
			((LD$BP_B - LD$BP_A) >= bins[n])),]$R2, na.rm=T),2)
		}

	}
	rm(LD)
}

# Average across chromosomes
data[,ncol(data)+1] <- round(rowMeans(data, na.rm=T),2)
rownames(data) <- bins

# Write to file
fname <- paste("../r2-pop-", pop, sep="")
write.table(data, fname, sep="\t", col.names=T, row.names=T, append=F, quote=F)



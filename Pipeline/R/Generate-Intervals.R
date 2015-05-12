setwd("/home/dave/Copy/HoneyBee/Analyses")

# Chromosome details
#@SQ     SN:1    LN:29893408     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa  M5:700cc70320e81afaad43cc7d0aeb58f6
#@SQ     SN:2    LN:15549267     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:001e7243bd8fee1e7e35573ab4cd026d
#@SQ     SN:3    LN:13234341     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:4cc93b0140c899cfdeb1f465bd5c73e7
#@SQ     SN:4    LN:12718334     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:539b1e7ded824b68ce978307d1e3eee0
#@SQ     SN:5    LN:14363272     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:bf9ae9e4a445d84edcef465c4e93b76a
#@SQ     SN:6    LN:18472937     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:070ca29f8c7497cda47a41c5149a2d72
#@SQ     SN:7    LN:13219345     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:136a93c49f7ff6f3929102a99c4acc53
#@SQ     SN:8    LN:13546544     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:9d367437d39e86f466e5108093cf17f8
#@SQ     SN:9    LN:11120453     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:104c9d42f3b16c2670d26424d26fc9f8
#@SQ     SN:10   LN:12965953     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:9fc075d9a331afa9f060a1f65360687a
#@SQ     SN:11   LN:14726556     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:334b137cbf4656e3f3bd20c4db016235
#@SQ     SN:12   LN:11902654     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:e8b161be17b7ad3705a8d43a4b779e17
#@SQ     SN:13   LN:10288499     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:ddb586b8eb1cb091bf523a2ba3fe89a9
#@SQ     SN:14   LN:10253655     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:0f212ad03266c452647bb4b0d4ad3d4e
#@SQ     SN:15   LN:10167229     UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:b36c5b62d77b704fec3cc9cffd0ffa96
#@SQ     SN:16   LN:7207165	UR:file:/home/dwragg/work/Analysis/GC/JFM10-TAGCTT-L001.fa	M5:24ccfd9be45ed26117bbc5c23e3a7d08

contig <- seq(1,16)
length <- c(29893408, 15549267, 13234341, 12718334, 14363272, 18472937, 13219345,
            13546544, 11120453, 12965953, 14726556, 11902654, 10288499, 10253655,
            10167229, 7207165)
contigs <- data.frame(contig, length, stringsAsFactors=F)

win_size <- 5000
win_step <- 1000

bed <- NULL
for(i in 1:nrow(contigs)){
  start <- as.integer(seq(1, contigs[i,2], by=win_step))
  chr <- rep(contigs[i,1], length(start))
  stop <- as.integer(start + win_size)
  stop[which(stop>contigs[i,2])] <- contigs[i,2]
  bed.tmp <- data.frame(chr, start, stop, stringsAsFactors=F)
  bed.tmp[nrow(bed.tmp),ncol(bed.tmp)] <- contigs[i,2]
  bed <- rbind(bed, bed.tmp)  
}
bed[,4] <- paste(bed[,1], ":", bed[,2], "-", bed[,3], sep="")

write.table(bed[,4], "intervals.list", sep="\t", row.names=F, col.names=F, append=F, quote=F)




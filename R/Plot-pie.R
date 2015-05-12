library(RPMG)

setwd("/home/dave/Copy/HoneyBee/Analyses/selection")

pdfname <- paste("Results/pie.png", sep="")
png(pdfname, height=175, width=175, units="mm", res=300)

bkgr="#FFFFFF"
frgr="#000000"

par(bg=bkgr, fg=frgr, col.lab=frgr, cex=0.8, col.axis=frgr,lwd=1)
vrow <- 2
vcol <- 2
par(mfrow=c(vrow,vcol), oma = c(1,1,1,1), mar = c(1,1,1,1))

cols <- pastel.colors(12, seed=13)

honey.slices <- c(1, 1, 1, 1, 1, 7, 8, 12, 20, 26)
honey.labels <- c("apoptotic process", "metabolism", "methylation", "carbohydrate metabolism", "biosynthesis", "mannose metabolism",
                  "chromatin modification", "ionotropic glutamate receptor signaling pathway","cobalamin transport",
                  "peptide cross-linking")
jelly.slices <- c(1, 1, 1, 1, 1, 2, 7, 13, 14, 16, 18)
jelly.labels <- c("apoptotic process", "cell adhesion", "metabolism",	"dephosphorylation",	"multicellular organismal development",
                  "biosynthesis", "protein homooligomerization", "oxygen transport", "steroid hormone mediated signaling pathway",
                  "carbohydrate metabolism", "peptidyl-lysine modification to peptidyl-hypusine")

pie(honey.slices, col=cols, labels=round(honey.slices/sum(honey.slices),3), cex=0.6)
title("Honey")
pie(jelly.slices, col=cols, labels=round(jelly.slices/sum(jelly.slices),3), cex=0.6)
title("Royal Jelly")

plot(seq(1,10,1), type="n", xaxt="n", yaxt="n", bty="n")
legend("top", bty="n", legend=rev(paste(honey.labels, " (", round(honey.slices/sum(honey.slices),3),")", sep="")), 
       cex=0.8, ncol=1, fill=rev(cols[1:length(honey.slices)]))

plot(seq(1,10,1), type="n", xaxt="n", yaxt="n", bty="n")
legend("top", bty="n", legend=rev(paste(jelly.labels, " (", round(jelly.slices/sum(jelly.slices),3),")", sep="")), 
       cex=0.8, ncol=1, fill=rev(cols[1:length(jelly.slices)]))



dev.off()





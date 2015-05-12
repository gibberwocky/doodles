# Currently the required GO package is not available for R 3.2.0, thus this script will not work

# source("http://bioconductor.org/biocLite.R")
# sudo apt-get install libxml2-dev
# biocLite("GSEABase")
# biocLite("GO")
library(GSEABase)

Dmel <- system.file("extdata", "/home/dave/Copy/HoneyBee/Annotation/Dmel-go-basic.obo", package="GSEABase")
GO <- getOBOCollection(Dmel)


GOTs <- c("GO:0007155", "GO:0007156", "GO:0007626", "GO:0008152", "GO:0015671", "GO:0034508", "GO:0050790", "GO:0043547",
           "GO:0044707", "GO:0016567", "GO:0032446", "GO:0042811", "GO:0006629", "GO:0007218", "GO:0006813", "GO:0006281",
           "GO:0051056", "GO:0006508", "GO:0055114", "GO:0042384", "GO:0042073", "GO:0006810", "GO:0016051", "GO:0070647",
           "GO:0007165", "GO:0006355", "GO:0044267", "GO:0008344", "GO:0006811", "GO:0006464", "GO:0036211", "GO:0006468",
           "GO:0030534", "GO:0055085", "GO:0034220")


tmp <-GOCollection(GOTs)
 


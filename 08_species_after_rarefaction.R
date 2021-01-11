# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


OTUs <- read.table("data/01_01_rarefied_otus_ed.txt", h=T)

library(vegan)

specnumber <- specnumber(OTUs[,6:4623], OTUs$treatment)
write.table(specnumber, "out/08_specnumber_vegan.txt")

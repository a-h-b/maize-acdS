# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


OTUs<-read.table("data/01_01_rarefied_otus_ed.txt", h=T)
Soil<-as.factor(OTUs[,3])
Plant<-as.factor(OTUs[,4])
Depth<-as.factor(OTUs[,5])
Name<-as.factor(OTUs[,2])
Treatment<-as.factor(OTUs[,1])

library(vegan)

all <- adonis(formula = OTUs[, 6:4623] ~ Soil * Plant * Depth, method = "bray") 

OTUs_loam <- subset (OTUs, OTUs$Soil == "L")
Soil_L <- as.factor(OTUs_loam[,3])
Genotype_L <- as.factor(OTUs_loam[,4])
Depth_L <- as.factor(OTUs_loam[,5])
Name_L <- as.factor(OTUs_loam[,2])

OTUs_sand <- subset (OTUs, OTUs$Soil == "S")
Soils_S <- as.factor(OTUs_sand[,3])
Genotype_S <- as.factor(OTUs_sand[,4])
Depths_S <- as.factor(OTUs_sand[,5])
Names_S <- as.factor(OTUs_sand[,2])

library(vegan)
loam <- adonis(formula = OTUs_loam[, 6:4623] ~ Genotype_L * Depth_L, method = "bray") 

sand <- adonis (formula = OTUs_sand[, 6:4623] ~ Genotype_S * Depths_S, method = "bray")

all$aov.tab
test <- rbind(all$aov.tab, loam$aov.tab, sand$aov.tab)
write.table(test, "out/07_permanova_all_loam_sand.txt")

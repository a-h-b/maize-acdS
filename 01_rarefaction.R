# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


OTUs <- read.table("data/01_otu_ed.txt", h = T)
head(OTUs)
str(OTUs)
OTUs <- OTUs [,1:59]
str(OTUs)

OTUs <- OTUs[, -23] # L_WT_C_2

write.table(
  OTUs,
  "data/01_00_otus_unrarefactioned_tax.txt",
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE
)

tOTUs <- as.data.frame(t(OTUs))
Treatments <- row.names(tOTUs)
tOTUs <- cbind(Treatments, tOTUs)
write.table(
  tOTUs,
  "data/01_00_otus_transponiert_unrare.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)
str(tOTUs)
head(tOTUs)

tOTUs <- read.table("data/01_00_otus_transponiert_unrare.txt", h = T)
OTUs <- tOTUs
rm(tOTUs)

#check reads for rarefaction: 32242 L_WT_B_2
reads <- margin.table(as.matrix(OTUs[, 2:4619]), 1)
min(reads)
sort(reads, dec = TRUE)

# Rarefaction -------------------------------------------------------------
library(vegan)

OTUs_rarefied <- rrarefy(OTUs[, 2:4619], 46346)

write.table(
  OTUs_rarefied,
  "data/01_01_rarefied_otus.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

OTU_rare <- read.table("data/01_01_rarefied_otus.txt", h = T)
depth <- rep("D", 58)
OTU_rare <- cbind(Treatment = OTUs$Treatments, depth, OTU_rare)
str(OTU_rare)
head(OTU_rare)
library(dplyr)
library(tidyverse)
library(tidyr)
OTU_rare_ed <- OTU_rare %>% 
  separate(col = Treatment, into = c("Soil", "Plant", "Letter","Depthr"), sep = "_") %>% 
  unite(Depth, depth, Depthr, sep = "", remove = FALSE) 

treatment <- OTUs$Treatments
OTU_rare_ed <- cbind(treatment, OTU_rare_ed)
OTU_rare_ed <- OTU_rare_ed[,-c(4,6,7)]
OTU_rare_ed <- OTU_rare_ed %>% 
  unite(Name, Soil:Depth, sep = "_", remove = FALSE)
Name <- OTU_rare_ed$Name

write.table(
  OTU_rare_ed,
  "data/01_01_rarefied_otus_ed.txt",
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE
)

# rarecurve

OTUs <- cbind(Name, OTUs[, 2:4619])
Treatment <- as.factor(OTUs[, 1])

str(OTUs)

# Soil Plant --------------------------------------------------------------
colour = rep(NA, length = length(OTUs$Name))
colour[which(OTUs$Name == "L_RTH_D1")] = "red3"
colour[which(OTUs$Name == "L_RTH_D2")] = "red3"
colour[which(OTUs$Name == "L_RTH_D3")] = "red3"
colour[which(OTUs$Name == "L_WT_D1")] = "blue3"
colour[which(OTUs$Name == "L_WT_D2")] = "blue3"
colour[which(OTUs$Name == "L_WT_D3")] = "blue3"
colour[which(OTUs$Name == "S_RTH_D1")] = "darkorange2"
colour[which(OTUs$Name == "S_RTH_D2")] = "darkorange2"
colour[which(OTUs$Name == "S_RTH_D3")] = "darkorange2"
colour[which(OTUs$Name == "S_WT_D1")] = "deepskyblue"
colour[which(OTUs$Name == "S_WT_D2")] = "deepskyblue"
colour[which(OTUs$Name == "S_WT_D3")] = "deepskyblue"

tiff("out/01_rarecurve_1.tiff", units="in", width=6, height=5, res=300)
rarecurve(OTUs[,2:4619], step=100, label=FALSE, col=colour, main="", xlab = "Reads", ylab = "OTUs", cex.lab = 1.2)     #Package VEGAN
legend(legend=c("L WT", "L RTH", 
                "S WT", "S RTH"),
       "bottomright", 
       bty="n",
       col=c("red3", "darkorange2", "blue3", "deepskyblue"),
       lty=1,
       pt.cex=1.5,
       cex=0.85,
       ncol=2)
arrows(x0=22912, y0=-50, y1=600, angle=90, length=0)
dev.off()


# Soil Depth --------------------------------------------------------------
colour = rep(NA, length = length(OTUs$Name))
colour[which(OTUs$Name == "L_RTH_D1")] = "orange"
colour[which(OTUs$Name == "L_RTH_D2")] = "red"
colour[which(OTUs$Name == "L_RTH_D3")] = "darkred"
colour[which(OTUs$Name == "L_WT_D1")] = "orange"
colour[which(OTUs$Name == "L_WT_D2")] = "red"
colour[which(OTUs$Name == "L_WT_D3")] = "darkred"
colour[which(OTUs$Name == "S_RTH_D1")] = "skyblue2"
colour[which(OTUs$Name == "S_RTH_D2")] = "blue2"
colour[which(OTUs$Name == "S_RTH_D3")] = "navy"
colour[which(OTUs$Name == "S_WT_D1")] = "skyblue2"
colour[which(OTUs$Name == "S_WT_D2")] = "blue"
colour[which(OTUs$Name == "S_WT_D3")] = "navy"

tiff("out/01_rarecurve_2.tiff", units="in", width=6, height=5, res=300)
rarecurve(OTUs[,2:4619], step=100, label=FALSE, col=colour, main="", xlab = "Reads", ylab = "OTUs", cex.lab = 1.2)     #Package VEGAN
legend(legend=c("L Depth 1", "L Depth 2", "L Depth 3",
                "S Depth 1", "S Depth 2", "S Depth 3"),
       "bottomright", 
       bty="n",
       col=c("orange", "red", "darkred", "skyblue2", "blue", "navy"),
       pch=15,
       pt.cex=1.5,
       cex=0.75,
       ncol=2)
arrows(x0=22912, y0=-50, y1=600, angle=90, length=0)
dev.off()


# All ---------------------------------------------------------------------
colour = rep(NA, length = length(OTUs$Name))
colour[which(OTUs$Name == "L_RTH_D1")] = "red3"
colour[which(OTUs$Name == "L_RTH_D2")] = "red3"
colour[which(OTUs$Name == "L_RTH_D3")] = "red3"
colour[which(OTUs$Name == "L_WT_D1")] = "blue3"
colour[which(OTUs$Name == "L_WT_D2")] = "blue3"
colour[which(OTUs$Name == "L_WT_D3")] = "blue3"
colour[which(OTUs$Name == "S_RTH_D1")] = "darkorange2"
colour[which(OTUs$Name == "S_RTH_D2")] = "darkorange2"
colour[which(OTUs$Name == "S_RTH_D3")] = "darkorange2"
colour[which(OTUs$Name == "S_WT_D1")] = "deepskyblue"
colour[which(OTUs$Name == "S_WT_D2")] = "deepskyblue"
colour[which(OTUs$Name == "S_WT_D3")] = "deepskyblue"

tiff("out/01_rarecurve_3.tiff", units="in", width=6, height=5, res=300)
rarecurve(OTUs[,2:4619], step=100, label=FALSE, col=colour, main="", xlab = "Reads", ylab = "OTUs",lty=c(1,2,3), cex.lab = 1.2)     #Package VEGAN
legend(legend=c("L WT D1", "L WT D2", "L WT D3", "L RTH D1", "L RTH D2", "L RTH D3", 
                "S WT D1", "S WT D2", "S WT D3", "S RTH D1", "S RTH D2", "S RTH D3"),
       "bottomright", 
       bty = "n",
       col=c(rep("blue",3), rep("red",3), rep("deepskyblue",3), rep("darkorange2",3)),
       lty=c(1,2,3),
       #lty.cex=1.5,
       cex=0.85,
       ncol=2)
arrows(x0=22912, y0=-50, y1=600, angle=90, length=0)
dev.off()


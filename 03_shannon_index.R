# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


OTUs <- read.table("data/01_01_rarefied_otus_ed.txt", h = T)
Precol <- read.table("data/01_01_rarefied_otus_ed.txt")[, 1:5]
str(OTUs)
Soil <- as.factor(OTUs[, 3])
Plant <- as.factor(OTUs[, 4])
Depth <- as.factor(OTUs[, 5])
Name <- as.factor(OTUs[, 2])
Treatment <- as.factor(OTUs[, 1])

#Shannon-Index: Package:VEGAN

library(vegan)

SI <- diversity(OTUs[, 6:4623], index = "shannon")

SI_table <- cbind(Precol, SI)

write.table(
  SI_table,
  "data/03_shannon_index.txt",
  sep = "\t",
  row.name = FALSE,
  col.name = TRUE
)

#bei Verwendung im Excel: im .txt nicht vergessen . durch , zu ersetzen :)

str(SI_table)
SI <- SI_table

Soil <- as.factor(SI[, 3])
Plant <- as.factor(SI[, 4])
Depth <- as.factor(SI[, 5])
Name <- as.factor(SI[, 2])
Treatment <- as.factor(SI[, 1])

L_WT_1 <- SI[grep("L_WT.*_1", SI$treatment), ]
L_WT_1_SI <- mean(L_WT_1[, 6])
L_WT_1_SD <- sd(L_WT_1[, 6])

L_WT_2 <- SI[grep("L_WT.*_2", SI$treatment), ]
L_WT_2_SI <- mean(L_WT_2[, 6])
L_WT_2_SD <- sd(L_WT_2[, 6])

L_WT_3 <- SI[grep("L_WT.*_3", SI$treatment), ]
L_WT_3_SI <- mean(L_WT_3[, 6])
L_WT_3_SD <- sd(L_WT_3[, 6])

L_RTH_1 <- SI[grep("L_RTH.*_1", SI$treatment), ]
L_RTH_1_SI <- mean(L_RTH_1[, 6])
L_RTH_1_SD <- sd(L_RTH_1[, 6])

L_RTH_2 <- SI[grep("L_RTH.*_2", SI$treatment), ]
L_RTH_2_SI <- mean(L_RTH_2[, 6])
L_RTH_2_SD <- sd(L_RTH_2[, 6])

L_RTH_3 <- SI[grep("L_RTH.*_3", SI$treatment), ]
L_RTH_3_SI <- mean(L_RTH_3[, 6])
L_RTH_3_SD <- sd(L_RTH_3[, 6])

S_WT_1 <- SI[grep("S_WT.*_1", SI$treatment), ]
S_WT_1_SI <- mean(S_WT_1[, 6])
S_WT_1_SD <- sd(S_WT_1[, 6])

S_WT_2 <- SI[grep("S_WT.*_2", SI$treatment), ]
S_WT_2_SI <- mean(S_WT_2[, 6])
S_WT_2_SD <- sd(S_WT_2[, 6])

S_WT_3 <- SI[grep("S_WT.*_3", SI$treatment), ]
S_WT_3_SI <- mean(S_WT_3[, 6])
S_WT_3_SD <- sd(S_WT_3[, 6])

S_RTH_1 <- SI[grep("S_RTH.*_1", SI$treatment), ]
S_RTH_1_SI <- mean(S_RTH_1[, 6])
S_RTH_1_SD <- sd(S_RTH_1[, 6])

S_RTH_2 <- SI[grep("S_RTH.*_2", SI$treatment), ]
S_RTH_2_SI <- mean(S_RTH_2[, 6])
S_RTH_2_SD <- sd(S_RTH_2[, 6])

S_RTH_3 <- SI[grep("S_RTH.*_3", SI$treatment), ]
S_RTH_3_SI <- mean(S_RTH_3[, 6])
S_RTH_3_SD <- sd(S_RTH_3[, 6])
#---

#library(PMCMRplus)
library(agricolae)

#für boxplot: 
shapiro.test(SI_table[1:29, 6]) #Loam normal verteilt
shapiro.test(SI_table[30:58, 6]) #Sand nicht normal verteilt

SI_loam <- SI[1:29, ]
Soil_loam <- as.factor(SI_loam[, 3])
Plant_loam <- as.factor(SI_loam[, 4])
Depth_loam <- as.factor(SI_loam[, 5])
Name_loam <- as.factor(SI_loam[, 2])
Treatment_loam <- as.factor(SI_loam[, 1])

aov_loam <- aov(SI_loam[, 6] ~ Name_loam, SI_loam)
summary(aov_loam) #0.000238 ***
aov_loam
groups_loam <- HSD.test(aov_loam, "Name_loam")
groups_loam$groups

SI_sand <- SI[30:58, ]
Soil_sand <- as.factor(SI_sand[, 3])
Plant_sand <- as.factor(SI_sand[, 4])
Depth_sand <- as.factor(SI_sand[, 5])
Name_sand <- as.factor(SI_sand[, 2])
Treatment_sand <- as.factor(SI_sand[, 1])
krusk_sand <-
  kruskal(SI_sand[, 6],
          Name_sand,
          p.adj = c("BH"),
          group = TRUE)
summary(krusk_sand)
krusk_sand$groups
krusk_sand$statistics
comp <- krusk_sand$comparison
write.table(comp, "out/20201216_Shannon-Index_pvalues_comparison.txt", row.names = TRUE, col.names = TRUE)

groups_loam

#Boxplot both

tiff("out/03_SI_loam_sand.tiff", units="in", width=12, height=6, res=300)
par(oma = c(0,0,0,0), mfrow = c(1,2))
loam <- boxplot(L_WT_1[,6], L_WT_2[,6], L_WT_3[,6], L_RTH_1[,6], L_RTH_2[,6], L_RTH_3[,6], main="", 
                ylab="", xlab="", ylim=c(4.5,6.3),
                col=c("blue3","blue3","blue3","red3","red3","red3"), 
                names=c("D1","D2","D3","D1", "D2","D3"), cex.lab = 1.5)
legend("topright", legend=c("WT", "RTH"), title="Plants:", fill=c("blue3","red3"), cex=1, bty='n')
mtext("Shannon-index", side=2, line=2.5, cex=1.5)
text(0, 6.3, "a)", xpd = T, cex = 1.5, font = 2)
text(c(1, 2, 3, 4, 5, 6), loam$stats[5,]+0.08, labels=c("b", "a", "ab", "a", "a", "ab"), cex = 1.5)
sand <- boxplot(S_WT_1[,6], S_WT_2[,6], S_WT_3[,6], S_RTH_1[,6], S_RTH_2[,6], S_RTH_3[,6], main="", 
                ylab="", xlab="", ylim=c(4.5,6.3),col=c("deepskyblue","deepskyblue","deepskyblue","darkorange2","darkorange2","darkorange2"), 
                names=c("D1","D2","D3","D1", "D2","D3"), cex.lab = 1.5)
legend("topright", legend=c("WT", "RTH"), title="Plants:", fill=c("deepskyblue","darkorange2"), cex=1, bty='n')
text(c(1, 2, 3, 4, 5, 6), sand$stats[5,]+0.08, labels=c("ab", "ab", "a", "b", "ab", "ab"), cex = 1.5)
text(0, 6.3, "b)", xpd = T, cex = 1.5, font = 2)
dev.off()
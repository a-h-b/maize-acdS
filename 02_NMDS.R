# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


OTUs <- read.table("data/01_01_rarefied_otus_ed.txt", h = T)


Soil <- as.factor(OTUs[, 3])
Plant <- as.factor(OTUs[, 4])
Depth <- as.factor(OTUs[, 5])
Name <- as.factor(OTUs[, 2])
Treatment <- as.factor(OTUs[, 1])

library(vegan)
tiff("out/02_stressplot_all.tiff", units="in", width=6, height=6, res=300)
NMDS_all <- metaMDS(OTUs[, 6:4623], k = 2)
stressplot_all <- stressplot(NMDS_all)
NMDS_all$stress
mtext("stress = 0,11")
dev.off()


# D1 ----------------------------------------------------------------------


OTUs_D1 <- OTUs[which((regexpr(".*1.*", OTUs$Depth)) == 1), ]
write.table(OTUs_D1, "13_otus_depth_1.txt", sep = "\t")
OTUs_D1 <- read.table("13_otus_depth_1.txt", h = T)
Soil_D1 <- as.factor(OTUs_D1[, 1])
Plant_D1 <- as.factor(OTUs_D1[, 2])
Depth_D1 <- as.factor(OTUs_D1[, 3])
Name_D1 <- as.factor(OTUs_D1[, 4])
Treatment_D1 <- as.factor(OTUs_D1[, 5])

tiff("13_stressplot_D1.tiff", units="in", width=6, height=6, res=300)
NMDS_D1 <- metaMDS(OTUs_D1[, 6:4088], k = 2)
stressplot_D1 <- stressplot(NMDS_D1)
NMDS_D1$stress
mtext("stress = 0,06")
dev.off()


# D2 ----------------------------------------------------------------------


OTUs_D2 <- OTUs[which((regexpr(".*2.*", OTUs$Depth)) == 1), ]
write.table(OTUs_D2, "13_otus_depth_2.txt", sep = "\t")
OTUs_D2 <- read.table("13_otus_depth_2.txt", h = T)
Soil_D2 <- as.factor(OTUs_D2[, 1])
Plant_D2 <- as.factor(OTUs_D2[, 2])
Depth_D2 <- as.factor(OTUs_D2[, 3])
Name_D2 <- as.factor(OTUs_D2[, 4])
Treatment_D2 <- as.factor(OTUs_D2[, 5])

tiff("13_stressplot_D2.tiff", units="in", width=6, height=6, res=300)
NMDS_D2 <- metaMDS(OTUs_D2[, 6:4088], k = 2)
stressplot_D2 <- stressplot(NMDS_D2)
NMDS_D2$stress
mtext("stress = 0,11")
dev.off()


# D3 ----------------------------------------------------------------------


OTUs_D3 <- OTUs[which((regexpr(".*3.*", OTUs$Depth)) == 1), ]
write.table(OTUs_D3, "13_otus_depth_3.txt", sep = "\t")
OTUs_D3 <- read.table("13_otus_depth_3.txt", h = T)
Soil_D3 <- as.factor(OTUs_D3[, 1])
Plant_D3 <- as.factor(OTUs_D3[, 2])
Depth_D3 <- as.factor(OTUs_D3[, 3])
Name_D3 <- as.factor(OTUs_D3[, 4])
Treatment_D3 <- as.factor(OTUs_D3[, 5])

tiff("13_stressplot_D3.tiff", units="in", width=6, height=6, res=300)
NMDS_D3 <- metaMDS(OTUs_D3[, 6:4088], k = 2)
stressplot_D3 <- stressplot(NMDS_D3)
NMDS_D3$stress
mtext("stress=0,07")
dev.off()


# WT ----------------------------------------------------------------------


OTUs_WT <- OTUs[which((regexpr(".*WT.*", OTUs$Plant)) == 1), ]
write.table(OTUs_WT, "13_otus_WT.txt", sep = "\t")
OTUs_WT <- read.table("13_otus_WT.txt", h = T)
Soil_WT <- as.factor(OTUs_WT[, 1])
Plant_WT <- as.factor(OTUs_WT[, 2])
Depth_WT <- as.factor(OTUs_WT[, 3])
Name_WT <- as.factor(OTUs_WT[, 4])
Treatment_WT <- as.factor(OTUs_WT[, 5])

tiff("13_stressplot_WT.tiff", units="in", width=6, height=6, res=300)
NMDS_WT <- metaMDS(OTUs_WT[, 6:4088], k = 2)
stressplot_WT <- stressplot(NMDS_WT)
NMDS_WT$stress
mtext("stress = 0,09")
dev.off()

# RTH ---------------------------------------------------------------------

OTUs_RTH <- OTUs[which((regexpr(".*RTH.*", OTUs$Plant)) == 1), ]
write.table(OTUs_RTH, "13_otus_RTH.txt", sep = "\t")
OTUs_RTH <- read.table("13_otus_RTH.txt", h = T)
Soil_RTH <- as.factor(OTUs_RTH[, 1])
Plant_RTH <- as.factor(OTUs_RTH[, 2])
Depth_RTH <- as.factor(OTUs_RTH[, 3])
Name_RTH <- as.factor(OTUs_RTH[, 4])
Treatment_RTH <- as.factor(OTUs_RTH[, 5])



tiff("13_stressplot_RTH.tiff", units="in", width=6, height=6, res=300)
NMDS_RTH <- metaMDS(OTUs_RTH[, 6:4088], k = 2)
stressplot_RTH <- stressplot(NMDS_RTH)
NMDS_RTH$stress
mtext("stress = 0,08")
dev.off()


# Sand --------------------------------------------------------------------

OTUs_Sand <- OTUs[which((regexpr("S.*", OTUs$Soil)) == 1), ]
write.table(OTUs_Sand, "data/02_otus_sand_rare.txt", sep = "\t")
OTUs_Sand <- read.table("data/02_otus_sand_rare.txt", h = T)
Soil_Sand <- as.factor(OTUs_Sand[, 3])
Plant_Sand <- as.factor(OTUs_Sand[, 4])
Depth_Sand <- as.factor(OTUs_Sand[, 5])
Name_Sand <- as.factor(OTUs_Sand[, 2])
Treatment_Sand <- as.factor(OTUs_Sand[, 1])



tiff("out/02_stressplot_S.tiff", units="in", width=6, height=6, res=300)
NMDS_Sand <- metaMDS(OTUs_Sand[, 6:4623], k = 2)
stressplot_Sand <- stressplot(NMDS_Sand)
NMDS_Sand$stress
mtext("stress = 0,13")
dev.off()

# Loam --------------------------------------------------------------------
OTUs_Loam <- OTUs[which((regexpr("L.*", OTUs$Soil)) == 1), ]
write.table(OTUs_Loam, "data/02_otus_loam_rare.txt", sep = "\t")
OTUs_Loam <- read.table("data/02_otus_loam_rare.txt", h = T)
Soil_Loam <- as.factor(OTUs_Loam[, 3])
Plant_Loam <- as.factor(OTUs_Loam[, 4])
Depth_Loam <- as.factor(OTUs_Loam[, 5])
Name_Loam <- as.factor(OTUs_Loam[, 2])
Treatment_Loam <- as.factor(OTUs_Loam[, 1])

tiff("out/02_stressplot_L.tiff", units="in", width=6, height=6, res=300)
NMDS_Loam <- metaMDS(OTUs_Loam[, 6:4623], k = 2)
stressplot_Loam <- stressplot(NMDS_Loam)
NMDS_Loam$stress
mtext("stress = 0,12")
dev.off()

# NMDS GGPLOT2 ------------------------------------------------------------
library(ggplot2)

# Loam --------------------------------------------------------------------

library(ggplot2)

round(NMDS_Loam$stress, 2)

tiff("out/02_NMDS_L.tiff", units="in", width=6, height=5, res=300)
scrs_l <- scores(NMDS_Loam, display = 'sites')
scrs_l <- cbind(as.data.frame(scrs_l), treatment = OTUs_Loam$Name)
cent_l <- aggregate(cbind(NMDS1, NMDS2) ~ treatment, data = scrs_l, FUN = mean)
segs_l <- merge(scrs_l, setNames(cent_l, c('treatment', 'oNMDS1', 'oNMDS2')),
                by = 'treatment', sort = FALSE)
gr.use_l <- factor(scrs_l$treatment)
col.gr_l <- c("#74C476","#31A354","#006D2C","#74A9CF","#2B8CBE","#045A8D")
color_l <- c(2, 2, 3, 4, 5, 6, 7)
color_transparent_l <- adjustcolor(col.gr_l, alpha.f = 1)
color_transparent_l[gr.use_l]
mult_l <- 1
mult1_l <- 1
loam.ggplot <- ggplot(data = scrs_l, aes(y = NMDS2, x = NMDS1)) +
  #coord_cartesian(xlim = c(-1.2, 2)) +
  #coord_cartesian(ylim = c(-1.6, 1)) +
  geom_segment(data = segs_l, mapping = aes(xend = oNMDS1, yend = oNMDS2,
                                            colour = treatment), alpha = 0.5, linetype = 3) +
  geom_point(data = cent_l, size = 2, aes(colour = treatment, shape = treatment), alpha = 5/10) +
  geom_point(aes(color = scrs_l$treatment, shape = scrs_l$treatment), size = 5) +
  scale_color_manual(labels = c("RTH D1", "RTH D2", "RTH D3", "WT D1", "WT D2", "WT D3"), name="",
                     values = c("L_RTH_D1" = "red3", "L_RTH_D2" = "red3", "L_RTH_D3"="red3",  
                                "L_WT_D1" = "blue3", "L_WT_D2" = "blue3", "L_WT_D3" = "blue3")) +
  scale_shape_manual(labels = c("RTH D1", "RTH D2", "RTH D3", "WT D1", "WT D2", "WT D3"), name="",
                     values = c("L_RTH_D1" = 16, "L_RTH_D2" = 13, "L_RTH_D3"=1,  
                                "L_WT_D1" = 16, "L_WT_D2" = 13, "L_WT_D3" = 1)) +
  annotate("text", Inf, -Inf, hjust = 1.05, vjust = -0.4,
           label = round(NMDS_Loam$stress, 2) , size = 5) +
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill = NA, size=0.5,
                                    linetype="solid"),
        panel.ontop = FALSE,
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 14))
loam.ggplot
dev.off()

# Sand --------------------------------------------------------------------

round(NMDS_Sand$stress, 2)

tiff("out/02_NMDS_S.tiff", units="in", width=6, height=5, res=300)
scrs_s <- scores(NMDS_Sand, display = 'sites')
scrs_s <- cbind(as.data.frame(scrs_s), treatment = OTUs_Sand$Name)
cent_s <- aggregate(cbind(NMDS1, NMDS2) ~ treatment, data = scrs_s, FUN = mean)
segs_s <- merge(scrs_s, setNames(cent_s, c('treatment', 'oNMDS1', 'oNMDS2')),
                by = 'treatment', sort = FALSE)
gr.use_s <- factor(scrs_s$treatment)
col.gr_s <- c("#74C476","#31A354","#006D2C","#74A9CF","#2B8CBE","#045A8D")
color_s <- c(2, 2, 3, 4, 5, 6, 7)
color_transparent_s <- adjustcolor(col.gr_s, alpha.f = 1)
color_transparent_s[gr.use_s]
mult_s <- 1
mult1_s <- 1
sand.ggplot <- ggplot(data = scrs_s, aes(y = NMDS2, x = NMDS1)) +
  #coord_cartesian(xlim = c(-1.2, 2)) +
  #coord_cartesian(ylim = c(-1.6, 1)) +
  geom_segment(data = segs_s, mapping = aes(xend = oNMDS1, yend = oNMDS2,
                                            colour = treatment), alpha = 0.5, linetype = 3) +
  geom_point(data = cent_s, size = 2, aes(colour = treatment, shape = treatment), alpha = 5/10) +
  geom_point(aes(color = scrs_s$treatment, shape = scrs_s$treatment), size = 5) +
  scale_color_manual(labels = c("RTH D1", "RTH D2", "RTH D3", "WT D1", "WT D2", "WT D3"), name="",
                     values = c("S_RTH_D1" = "darkorange2", "S_RTH_D2" = "darkorange2", "S_RTH_D3"="darkorange2",  
                                "S_WT_D1" = "deepskyblue", "S_WT_D2" = "deepskyblue", "S_WT_D3" = "deepskyblue")) +
  scale_shape_manual(labels = c("RTH D1", "RTH D2", "RTH D3", "WT D1", "WT D2", "WT D3"), name="",
                     values = c("S_RTH_D1" = 16, "S_RTH_D2" = 13, "S_RTH_D3"=1,  
                                "S_WT_D1" = 16, "S_WT_D2" = 13, "S_WT_D3" = 1)) +
  annotate("text", Inf, -Inf, hjust = 1.05, vjust = -0.4,
           label = round(NMDS_Sand$stress, 2) , size = 5) +
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill = NA, size=0.5,
                                    linetype="solid"),
        panel.ontop = FALSE,
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.line = element_line(size = 0.5, colour = "black"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.text = element_text(size = 14))
sand.ggplot
dev.off()

# NMDS both ---------------------------------------------------------------

library(ggpubr)

tiff("out/02_NMDS_loam_sand.tiff", units="in", width=13, height=7, res=300)
ggarrange( loam.ggplot, sand.ggplot,
           labels = c("a)", "b)"),
           ncol = 2, nrow = 1, common.legend = FALSE, legend = "bottom")
dev.off()

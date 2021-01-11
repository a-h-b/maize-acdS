# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


pacman::p_load(betapart, conflicted, vegan)

comm <- read.table("data/01_01_rarefied_otus_ed.txt")
treat <- as.factor(comm[,2])

res_betaab <- beta.pair.abund(betapart.core.abund(comm[,6:4623]), index.family = "bray")
abu <- betadisper(res_betaab[[3]], treat)
abu
plot (abu)
scores_abu <- scores(abu, display = "sites")
eigens <- eigenvals(abu)
summary(abu, scaling = 1)
summary(abu, scaling = 2)
plot(abu)
site.scaling1 <- summary(scores_abu, scaling = 1)
env.scaling2 <- summary(scores_abu, scaling = 2)
eigenirg <- eigens/sum(eigens)
round(eigenirg[1]*100, 2)
round(eigenirg[2]*100, 2)
abu$group

tiff("out/10_PCA_BD.tiff", units="in", width=7, height=5, res=300)
plot(abu,
     main = "", 
     xlab = "", 
     ylab = "",
     ellipse = FALSE, sub = "",  label = FALSE, cex = 1.5, hull = FALSE,
     col = c(rep("red3",3), rep("blue3",3), rep("darkorange2", 3), rep("deepskyblue", 3)), 
     pch = c(rep(c(16,13,1),7),16,1, rep(c(16,13,1),6),16,13,rep(c(16,13,1),5)),
     seg.col = "black", seg.lty = 3, seg.lwd = 0.1)
mtext(paste("PCA 2 (", round(eigenirg[2]*100, 2), "%)", sep = ""), side = 2, line = 2.5)
mtext(paste("PCA 1 (", round(eigenirg[1]*100, 2), "%)", sep = ""), side = 1, line = 2.5)
arrows(-1,0,1, length = 1, angle = 90, lty = 2)
arrows(0,1,0,-1, length = 1, angle = 180, lty = 2)
legend(-0.25, 0.38, legend=c( "L WT 1", "L WT 2", "L WT 3","L RTH 1", "L RTH 2", "L RTH 3", "S WT 1", "S WT 2", "S WT 3", "S RTH 1", "S RTH 2", "S RTH 3"),
       bty = "n", bg = "white", xpd = T,
       col= c  ("blue3", "blue3", "blue3","red3", "red3", "red3", "deepskyblue", "deepskyblue", "deepskyblue", "darkorange2", "darkorange2", "darkorange2"),
       pch=(rep(c(16,13,1),4)),
       pt.cex=1.5,
       cex=0.75,
       ncol=4)
dev.off()

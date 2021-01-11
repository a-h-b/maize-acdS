# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


OTUs <- read.table("data/04_Genus_table_rare_sum.txt",h=T)[,c(1,6:79)]
precol<-read.table("data/04_Genus_table_rare_sum.txt",h=T)[,1:5]


# Vorgeplänkel ------------------------------------------------------------

L_WT_1<-OTUs[grep("L_WT.*_1", OTUs$treatment),]
L_WT_1_S<-colMeans(L_WT_1[,2:75])

L_WT_2<-OTUs[grep("L_WT.*_2", OTUs$treatment),]
L_WT_2_S<-colMeans(L_WT_2[,2:75])

L_WT_3<-OTUs[grep("L_WT.*_3", OTUs$treatment),]
L_WT_3_S<-colMeans(L_WT_3[,2:75])

L_RTH_1<-OTUs[grep("L_RTH.*_1", OTUs$treatment),]
L_RTH_1_S<-colMeans(L_RTH_1[,2:75])

L_RTH_2<-OTUs[grep("L_RTH.*_2", OTUs$treatment),]
L_RTH_2_S<-colMeans(L_RTH_2[,2:75])

L_RTH_3<-OTUs[grep("L_RTH.*_3", OTUs$treatment),]
L_RTH_3_S<-colMeans(L_RTH_3[,2:75])

S_WT_1<-OTUs[grep("S_WT.*_1", OTUs$treatment),]
S_WT_1_S<-colMeans(S_WT_1[,2:75])

S_WT_2<-OTUs[grep("S_WT.*_2", OTUs$treatment),]
S_WT_2_S<-colMeans(S_WT_2[,2:75])

S_WT_3<-OTUs[grep("S_WT.*_3", OTUs$treatment),]
S_WT_3_S<-colMeans(S_WT_3[,2:75])

S_RTH_1<-OTUs[grep("S_RTH.*_1", OTUs$treatment),]
S_RTH_1_S<-colMeans(S_RTH_1[,2:75])

S_RTH_2<-OTUs[grep("S_RTH.*_2", OTUs$treatment),]
S_RTH_2_S<-colMeans(S_RTH_2[,2:75])

S_RTH_3<-OTUs[grep("S_RTH.*_3", OTUs$treatment),]
S_RTH_3_S<-colMeans(S_RTH_3[,2:75])

genera<-rbind(L_WT_1_S, L_WT_2_S, L_WT_3_S, L_RTH_1_S, L_RTH_2_S, L_RTH_3_S, "", S_WT_1_S, S_WT_2_S, S_WT_3_S, S_RTH_1_S, S_RTH_2_S, S_RTH_3_S)
genera_t<-rbind(L_WT_1_S, L_WT_2_S, L_WT_3_S, L_RTH_1_S, L_RTH_2_S, L_RTH_3_S, S_WT_1_S, S_WT_2_S, S_WT_3_S, S_RTH_1_S, S_RTH_2_S, S_RTH_3_S)

write.table(genera_t, "data/06_genus_table_rare_sum_per_treatment.txt", sep="\t", col.names=TRUE, row.names=TRUE)

OTUs <- read.table("data/06_genus_table_rare_sum_per_treatment.txt", h=T)
OTUs <- as.matrix(OTUs)

reads_ges <- margin.table(OTUs, 2)
reads_ges
sort(reads_ges, dec = T)
OTUs_margin <- rbind(OTUs, reads_ges)
OTUs_margin <- t(OTUs_margin)
OTUs_margin_sort <- OTUs_margin[order(-OTUs_margin[,13]),]
OTUs_margin_sort <- t(OTUs_margin_sort)

OTUs_sort <- as.data.frame(OTUs_margin_sort)

rest <- rowSums(OTUs_sort[,21:74])
colnames(OTUs_sort)

str(OTUs_sort)

Streptomyces <- as.numeric(OTUs_sort[,1])
Variovorax <- as.numeric(OTUs_sort[,2])
Tetrasphaera <- as.numeric(OTUs_sort[,3])
Burkholderia <- as.numeric(OTUs_sort[,4])
Amycolatopsis <- as.numeric(OTUs_sort[,5])
Acidovorax <- as.numeric(OTUs_sort[,6])
Actinobacteria_uncl <- as.numeric(OTUs_sort[,7])
Microbacterium <- as.numeric(OTUs_sort[,8])
Methylibium <- as.numeric(OTUs_sort[,9])
Agromyces <- as.numeric(OTUs_sort[,10])
Brevibacterium <- as.numeric(OTUs_sort[,11])
Ralstonia <- as.numeric(OTUs_sort[,12])
Modestobacter <- as.numeric(OTUs_sort[,13])
Paraburkholderia <- as.numeric(OTUs_sort[,14])
Phycicoccus <- as.numeric(OTUs_sort[,15])
Bacteria_uncl <- as.numeric(OTUs_sort[,16])
Methylobacterium <- as.numeric(OTUs_sort[,17])
Marmoricola  <- as.numeric(OTUs_sort[,18])
Pseudomonas <- as.numeric(OTUs_sort[,19])
Arthrobacter <- as.numeric(OTUs_sort[,20])


genera_sort <- cbind(Streptomyces, Tetrasphaera, Microbacterium, Agromyces, Brevibacterium, Phycicoccus,
                     Arthrobacter, Amycolatopsis, Modestobacter, Marmoricola,Actinobacteria_uncl, Variovorax,
                     Burkholderia, Acidovorax, Methylibium, Ralstonia, Paraburkholderia, Methylobacterium, 
                     Pseudomonas,Bacteria_uncl, rest)

rownames(genera_sort)<- c("1", "2", "3", "1", "2", "3",  "1", "2", "3", "1", "2", "3")

genera_p<-100*prop.table(genera_sort, 1)
margin.table(genera_p,1)
genera_p <- genera_p[-13,]
genera_p <- genera_p[,-21]

write.table(genera_p, "out/06_proportional_table_genera.txt", sep="\t", row.names = TRUE, col.names = TRUE)

genera_plot<-t(genera_p)
leer <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
genera_plot <- cbind(genera_plot[,1:6], leer, genera_plot[,7:12])
colnames(genera_plot) <- c("", "", "", "", "", "", "", "", "", "", "", "", "")


# Plot --------------------------------------------------------------------

tiff("out/06_genus_diagramm.tiff", units="in", width=6, height=6, res=300)
par(oma=c(0,0,0,11))
bp <- barplot(genera_plot, space=0.25, cex.names = 0.9,
              main="", xlab="", ylab="", las=0, sub="", ylim = c(0,100),
              #col=c("#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF", "#3690C0",  "#0570B0", "#045A8D", "#034A74",
              #      "#022B43", "black", 
              #      "#FFFFCC", "#FFEDA0", "#FED976",  "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026"))
              col = c("#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF", "#3690C0",  "#0570B0", "#045A8D", "#034A74",    #actino
                      "#022B43", "#0d1a26",                                                                                    #actino
                      "#FFFFCC", "#FFEDA0", "#FED976",  "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "black"))   #proteo + uncl
                    
text(bp, -10,c("", "WT", "", "", "RTH", "", "", "", "WT", "", "", "RTH", ""), xpd = TRUE)
text(c(3.75, 12.5), c(105, 105), c("Loam", "Sand"), xpd = TRUE)
mtext("Percentage of reads", side = 2, line = 2.5)
text(bp, -3, c("D1", "D2", "D3", "D1", "D2", "D3", "", "D1", "D2", "D3", "D1", 
               "D2", "D3"), xpd = TRUE, cex = 0.8)
arrows(0.25, 102, 7.5, length=0, xpd=T)
arrows(9, 102, 16.25, length=0, xpd=T)
arrows(0.25, -6.5, 3.75, length=0, xpd=T)
arrows(4, -6.5, 7.5, length=0, xpd=T)
arrows(9, -6.5, 12.5, length=0, xpd=T)
arrows(12.75, -6.5, 16.25, length=0, xpd=T)
legend(19,115, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="", text.font = c(1,3,3,3,3,1,3,3,3),
       legend = c("uncl. Bacteria ***"),
       fill =   c("black"))
legend(19, 106, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Pseudomonadales", text.font = c(3,3,3,3,3,1,3,3,3),
       legend = c("Pseudomonas *"),
       fill =   c("#BD0026"))
legend(19, 97, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Rhizobiales", text.font = c(3,3,3,3,3,1,3,3,3),
       legend = c("Methylobacterium **"),
       fill =   c("#E31A1C"))
legend(19,87.5, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Burkholderiales", text.font = c(3,3,3,3,3,3,3,3,3),
       legend = c("Paraburkholderia *", "Ralstonia ***", "Methylibium ***", "Acidovorax ***", "Burkholderia ***", "Variovorax ***"),
       fill =   c( "#FC4E2A", "#FD8D3C", "#FEB24C", "#FED976", "#FFEDA0", "#FFFFCC"))
legend(19,62.5, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="", text.font = c(1,3,3,3,3,1,3,3,3,3,3),
       legend = c("uncl. Actinobacteria ***"),
       fill =   c("black"))
legend(19,53.5, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Propionibacterales", text.font = c(3,3,3,3,3,1,3,3,3,3,3),
       legend = c("Marmoricola *"),
       fill =   c("#022B43"))
legend(19,44, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Geodermatophilales", text.font = c(3,3,3,3,3,1,3,3,3,3,3),
       legend = c("Modestobacter **"),
       fill =   c("#034A74"))
legend(19,35.5, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Pseudonocardiales", text.font = c(3,3,3,3,3,1,3,3,3,3,3),
       legend = c( "Amycolatopsis ***"),
       fill =   c( "#045A8D"))
legend(19,26.5, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Micrococcales", text.font = c(3,3,3,3,3,3,3,3,3,3,3),
       legend = c( "Arthrobacter", "Phycicoccus", "Brevibacterium **", "Agromyces ***", "Microbacterium ***", 
                   "Tetrasphaera ***"),
       fill =   c(  "#0570B0","#3690C0", "#74A9CF", "#A6BDDB", "#D0D1E6",
                   "#ECE7F2"))
legend(19,-3, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Streptomycetales", text.font = c(3,3,3,3,3,1,3,3,3,3,3),
       legend = c( "Streptomyces ***"),
       fill =   c( "#FFF7FB"))
text(17.5,84, "Proteobacteria", srt=90, xpd=T, cex = 0.85)
text(17.5,23, "Actinobacteria", srt=90, xpd=T, cex = 0.85)
arrows(18.5, -12.5, 18.5, 56, length=0, xpd=T)
arrows(16.75, -12.5, 16.75, 56, length=0, xpd=T)
arrows(16.75, -12.5, 18.5, -12.5, length=0, xpd=T)
arrows(16.75, 56, 18.5, 56, length=0, xpd=T)
arrows(18.5, 57.5, 18.5, 108, length=0, xpd=T)
arrows(16.75, 57.5, 16.75, 108, length=0, xpd=T)
arrows(18.5, 57.5, 18.5, 108, length=0, xpd=T)
arrows(16.75, 57.5, 18.5, 57.5, length=0, xpd=T)
arrows(16.75, 108, 18.5, 108, length=0, xpd=T)
dev.off()


# statistics --------------------------------------------------------------

stat <- cbind (precol, OTUs$Streptomyces, OTUs$Tetrasphaera, OTUs$Microbacterium, OTUs$Agromyces, OTUs$Brevibacterium, OTUs$Phycicoccus, OTUs$Arthrobacter, 
               OTUs$Amycolatopsis, OTUs$Modestobacter, OTUs$Marmoricola, OTUs$Actinobacteria_uncl, OTUs$Variovorax, OTUs$Burkholderia, OTUs$Acidovorax, 
               OTUs$Methylibium, OTUs$Ralstonia, OTUs$Paraburkholderia, OTUs$Methylobacterium, OTUs$Pseudomonas, OTUs$Bacteria_uncl)

Soil <- as.factor(stat[,3])
Plant <- as.factor(stat[,4])
Depth <- as.factor(stat[,5])
Name <- as.factor(stat[,2])
Treat <- as.factor(stat[,1])

library(agricolae)

top_20_kruskal <- kruskal (stat[,6], Name, p.adj = "BH",alpha = 0.05, group = TRUE)
top_20_kruskal$groups

f=file("out/06_results_kruskal_wallis_names_BH.txt","w")
for (i in 6:25){
        c = eval(stat[,i])
        p=kruskal(stat[,i], Name, p.adj = "BH", group = T)
        capture.output(print(eval(colnames(stat)[i])),file=f,append=TRUE)
        capture.output(p$statistics,file=f,append=TRUE)
}
close(f)

f=file("out/06_groups_kruskal_wallis_names_BH_statistics.txt","w")
for (i in 6:25){
        c = eval(stat[,i])
        p=kruskal(stat[,i], Name, p.adj = "BH", group = T)
        capture.output(print(eval(colnames(stat)[i])),file=f,append=TRUE)
        capture.output(p$groups,file=f,append=TRUE)
}
close(f)


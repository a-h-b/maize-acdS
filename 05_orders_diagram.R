# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


OTUs<-read.table("data/04_Order_table_rare_sum.txt",h=T)[,c(1,6:25)]
precol<-read.table("data/04_Order_table_rare_sum.txt",h=T)[,1:5]


str(OTUs)
L_WT_1<-OTUs[grep("L_WT.*_1", OTUs$treatment),]
L_WT_1_S<-colSums(L_WT_1[,2:21])

L_WT_2<-OTUs[grep("L_WT.*_2", OTUs$treatment),]
L_WT_2_S<-colSums(L_WT_2[,2:21])

L_WT_3<-OTUs[grep("L_WT.*_3", OTUs$treatment),]
L_WT_3_S<-colSums(L_WT_3[,2:21])

L_RTH_1<-OTUs[grep("L_RTH.*_1", OTUs$treatment),]
L_RTH_1_S<-colSums(L_RTH_1[,2:21])

L_RTH_2<-OTUs[grep("L_RTH.*_2", OTUs$treatment),]
L_RTH_2_S<-colSums(L_RTH_2[,2:21])

L_RTH_3<-OTUs[grep("L_RTH.*_3", OTUs$treatment),]
L_RTH_3_S<-colSums(L_RTH_3[,2:21])

S_WT_1<-OTUs[grep("S_WT.*_1", OTUs$treatment),]
S_WT_1_S<-colSums(S_WT_1[,2:21])

S_WT_2<-OTUs[grep("S_WT.*_2", OTUs$treatment),]
S_WT_2_S<-colSums(S_WT_2[,2:21])

S_WT_3<-OTUs[grep("S_WT.*_3", OTUs$treatment),]
S_WT_3_S<-colSums(S_WT_3[,2:21])

S_RTH_1<-OTUs[grep("S_RTH.*_1", OTUs$treatment),]
S_RTH_1_S<-colSums(S_RTH_1[,2:21])

S_RTH_2<-OTUs[grep("S_RTH.*_2", OTUs$treatment),]
S_RTH_2_S<-colSums(S_RTH_2[,2:21])

S_RTH_3<-OTUs[grep("S_RTH.*_3", OTUs$treatment),]
S_RTH_3_S<-colSums(S_RTH_3[,2:21])

Orders <- rbind(L_WT_1_S, L_WT_2_S, L_WT_3_S, L_RTH_1_S, L_RTH_2_S, L_RTH_3_S, "", S_WT_1_S, S_WT_2_S, S_WT_3_S, S_RTH_1_S, S_RTH_2_S, S_RTH_3_S)
Orders_t <- rbind(L_WT_1_S, L_WT_2_S, L_WT_3_S, L_RTH_1_S, L_RTH_2_S, L_RTH_3_S, S_WT_1_S, S_WT_2_S, S_WT_3_S, S_RTH_1_S, S_RTH_2_S, S_RTH_3_S)
write.table(Orders_t, "data/05_order_table_rare_sum_per_treatment.txt", sep="\t", col.names=TRUE, row.names=TRUE)
write.table(Orders, "data/05_order_table_rare_sum_per_treatment_plot.txt", sep="\t", col.names=TRUE, row.names=TRUE)
Orders<-read.table("data/05_order_table_rare_sum_per_treatment_plot.txt",h=T)
as.data.frame(Orders)

as.data.frame(Orders)
#actino<-as.numeric(Orders[,2])
#sorda<-as.numeric(Orders[,7])
#beta<-as.numeric(Orders[,1])
#gamma<-as.numeric(Orders[,5])
#alpha<-as.numeric(Orders[,4])
#pro_u<-as.numeric(Orders[,6])
#bac_u<-as.numeric(Orders[,3])
sum(Orders[-7,1])

Orders_sort<-cbind(Orders$Streptomycetales, Orders$Micrococcales, Orders$Pseudonocardiales, Orders$Geodermatophilales, Orders$Propionibacteriales, Orders$Corynebacteriales,
                   Orders$Micromonosporales, Orders$Nakamurellales, Orders$Actinobacteria_uncl, Orders$Pseudomonadales, Orders$Enterobacterales, Orders$Gammaproteobacteria_uncl,
                   Orders$Burkholderiales, Orders$Rhizobiales, Orders$Rhodobacterales, Orders$Rhodospirillales, Orders$Alphaproteobacteria_uncl, Orders$Proteobacteria_uncl,
                   Orders$Bacteria_uncl, Orders$Sordariales)
row.names(Orders_sort)<- c("1", "2", "3", "1", "2", "3", "", "1", "2", "3", "1", "2", "3")
Orders_p<-100*prop.table(Orders_sort, 1)
str(Orders_sort)
#überprüfen -> 12x 1
margin.table(Orders_p, 1)

Orders_plot<-t(Orders_p)
colnames(Orders_plot) <- c("", "", "", "", "", "", "", "", "", "", "", "", "")


# statistics --------------------------------------------------------------

library(agricolae)
stat <- cbind(precol, OTUs[,2:21])
Name <- as.factor(stat[,2])
order_kruskal <- kruskal (stat[,6], Name, p.adj = "BH",alpha = 0.05, group = TRUE)
order_kruskal$groups
order_kruskal$statistics
order_kruskal$parameters

f=file("out/05_Results_kruskal_wallis_names_BH.txt","w")
for (i in 6:25){
  c = eval(stat[,i])
  p=kruskal(stat[,i], Name, p.adj = "BH", group = T)
  capture.output(print(eval(colnames(stat)[i])),file=f,append=TRUE)
  capture.output(p$statistics,file=f,append=TRUE)
}
close(f)

f=file("out/05_groups_kruskal_wallis_names_BH.txt","w")
for (i in 6:25){
  c = eval(stat[,i])
  p=kruskal(stat[,i], Name, p.adj = "BH", group = T)
  capture.output(print(eval(colnames(stat)[i])),file=f,append=TRUE)
  capture.output(p$groups,file=f,append=TRUE)
}
close(f)
# plot --------------------------------------------------------------------

tiff("out/05_order_diagramm.tiff", units="in", width=6, height=6, res=300)
par(oma=c(0,0,0,11))
bp <- barplot(Orders_plot, space=0.25, cex.names = 0.9,
              main="", xlab="", ylab="", las=0, sub="",
              col = c( "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF", "#3690C0",  "#0570B0", "#045A8D", "#034A74", "#022B43", 
                       "#FFFFCC", "#FFEDA0", "#FED976",  
                       "#FEB24C", 
                       "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026",
                       "#800026", "black", 
                       "#31a354"))
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
text(17.75,98, "Eukaryote", srt=90, xpd=T, cex = 0.85)
text(17.75,30, "Prokaryote", srt=90, xpd=T, cex = 0.85)
arrows(18.5, -12.5, 18.5, 86.5, length=0, xpd=T)
arrows(16.75, -12.5, 16.75, 86.5, length=0, xpd=T)
arrows(16.75, -12.5, 18.5, -12.5, length=0, xpd=T)
arrows(16.75, 86.5, 18.5, 86.5, length=0, xpd=T)
arrows(18.5, 88, 18.5, 108, length=0, xpd=T)
arrows(16.75, 88, 16.75, 108, length=0, xpd=T)
arrows(18.5, 88, 18.5, 108, length=0, xpd=T)
arrows(16.75, 88, 18.5, 88, length=0, xpd=T)
arrows(16.75, 108, 18.5, 108, length=0, xpd=T)
legend(19,108, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Sordariomycetes", text.font = c(3,3,3,3,3,1,3,3,3),
       legend = c("Sordariales"),
       fill =   c("#31a354"))
legend(19, 93.5, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title=" ", text.font = c(1,1,3,3,3,1,3,3,3),
       legend = c("uncl. Bacteria ***"),
       fill =   c("black"))
legend(19, 89.5, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title=" ", text.font = c(1,1,3,3,3,1,3,3,3),
       legend = c("uncl. Proteobacteria"),
       fill =   c("#800026"))
legend(19, 80, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Alphaproteobacteria", text.font = c(1,3,3,3,3,3,3,3,3),
       legend = c( "unclassified", "Rhodospirillales *", "Rhodobacterales **", "Rhizobiales **"),
       fill =   c( "#BD0026", "#E31A1C", "#FC4E2A", "#FD8D3C"))
legend(19, 58, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Betaproteobacteria", text.font = c(3,3,3,3,3,3,3,3,3,3,3),
       legend = c("Burkhodleriales ***"),
       fill =   c("#FEB24C"))
legend(19, 48, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Gammaproteobacteria", text.font = c(1,3,3,3,3,3,3,3,3,3,3),
       legend = c("unclassified", "Enterobacterales", "Pseudomonadales *"),
       fill =   c("#FED976", "#FFEDA0", "#FFFFCC"))
legend(19, 30, cex = 0.85, xpd = NA, title.adj = 0,bty = "n", title="Actinobacteria", text.font = c(1,3,3,3,3,3,3,3,3,3,3),
       legend = c( "unclassified ***", "Nakamurellales", "Micromonosporales", "Corynebacteriales **", "Propionibacteriales **", "Geodermatophilales **", 
                   "Pseudonocardiales ***", "Micrococcales ***", "Streptomycetales ***"),
       fill =   c( "#022B43", "#034A74", "#045A8D",  "#0570B0", "#3690C0", "#74A9CF", "#A6BDDB", "#D0D1E6", "#ECE7F2"))

dev.off()


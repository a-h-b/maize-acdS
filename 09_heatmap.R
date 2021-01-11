# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


data <- read.table("data/06_genus_table_rare_sum_per_treatment.txt", h=T)
data_final <- cbind(data$Streptomyces, data$Tetrasphaera, data$Microbacterium, data$Agromyces, data$Brevibacterium, data$Phycicoccus,
                    data$Arthrobacter, data$Amycolatopsis, data$Modestobacter, data$Marmoricola, data$Actinobacteria_uncl, data$Variovorax,
                    data$Burkholderia, data$Acidovorax, data$Methylibium, data$Ralstonia, data$Paraburkholderia, data$Methylobacterium, 
                    data$Pseudomonas, data$Bacteria_uncl)
colnames(data_final) <- c("Streptomyces", "Tetrasphaera", "Microbacterium", "Agromyces", "Brevibacterium", "Phycicoccus",
                          "Arthrobacter", "Amycolatopsis", "Modestobacter", "Marmoricola","Actinobacteria_uncl", "Variovorax",
                          "Burkholderia", "Acidovorax", "Methylibium", "Ralstonia", "Paraburkholderia", "Methylobacterium", 
                          "Pseudomonas","Bacteria_uncl")
row.names(data_final) <- row.names(data)


L_WT_1 <- c("white", "white", "black", "white", "black", "white", "white", "white", "white", "black", 
            "white", "black", "black", "black", "black", "black", "black", "white", "black", "white")
L_WT_2 <- c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", 
            "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")
L_WT_3 <- c("white", "white", "black", "white", "black", "black", "black", "black", "black", "black", 
            "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")
L_RTH_1 <- c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", 
             "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")
L_RTH_2 <- c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", 
             "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")
L_RTH_3 <- c("white", "white", "black", "white", "black", "black", "white", "black", "white", "black", 
             "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")
S_WT_1 <- c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", 
            "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")
S_WT_2 <- c("black", "black", "white", "black", "black", "black", "black", "black", "black", "black", 
            "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")
S_WT_3 <- c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", 
            "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")
S_RTH_1 <- c("black", "black", "white", "black", "white", "black", "black", "black", "black", "black", 
             "black", "black", "black", "black", "black", "black", "black", "white", "white", "white")
S_RTH_2 <- c("black", "black", "black", "black", "black", "white", "black", "black", "black", "black", 
             "black", "black", "black", "black", "black", "black", "black", "black", "white", "black")
S_RTH_3 <- c("black", "black", "black", "black", "white", "black", "black", "black", "black", "white", 
             "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")

numbercolor <- as.matrix(cbind(L_WT_1, L_WT_2, L_WT_3, L_RTH_1, L_RTH_2, L_RTH_3, S_WT_1, S_WT_2, S_WT_3, S_RTH_1, S_RTH_2, S_RTH_3))

library(vegan)
plot_data <- decostand(t(data_final),"standardize",1)
colors <-  colorRampPalette(c("navy", "white", "firebrick3"))(50) 

groups <- read.table("data/06_groups.txt", skipNul = TRUE, h=T)
groups <- rbind(groups[20,], groups[19,], groups[18,], groups[17,], groups[16,], 
                groups[15,], groups[14,], groups[13,], groups[12,], groups[11,], 
                groups[10,], groups[9,], groups[8,], groups[7,], groups[6,], 
                groups[5,], groups[4,], groups[3,], groups[2,], groups[1,])
#groups <- rbind(groups[4:6,],groups[1:3,], groups[10:12,], groups[7:9,])
#rownames(groups) <- groups[,1]
#groups <- groups [,-1]

library(pheatmap)

tiff("out/09__heatmap_generas_rarefied_data.tiff", units="in", width=9, height=6, res=300)
plot.new()
pheatmap(plot_data, cluster_rows = FALSE, cluster_cols = FALSE, cex=1, cellwidth = 28, cellheight = 14, 
         angle = 90, border_color = "white", fontsize_number = 8.5, col = colors, fontsize = 13,
         #scale = "column", 
         breaks = seq(from=-2,to=2,length.out = 51),
         legend = TRUE,
         display_numbers = as.matrix(groups),
         number_color = numbercolor,
         legend_breaks = -2:2, 
         legend_labels = c( "", "-2",  " 0",  " 2", "Scaled rel. abundance:"),
         labels_col = c("L WT D1", "L WT D2", "L WT D3","L RTH D1", "L RTH D2", "L RTH D3", "S WT D1", "S WT D2", "S WT D3", "S RTH D1", "S RTH D2", "S RTH D3"),
         labels_row = c("un. Bacteria", expression(italic("Pseudomonas")),
                        expression(italic("Methylobacterium")), expression(italic("Paraburkholderia")),
                        expression(italic("Ralstonia")), expression(italic("Methylibium")), 
                        expression(italic("Acidovorax")), expression(italic("Burkholderia")),  
                        expression(italic("Variovorax")), "un. Actinobacteria",
                        expression(italic("Marmoricola")), expression(italic( "Modestobacter")), 
                        expression(italic("Amycolatopsis")), expression(italic("Arthrobacter")), 
                        expression(italic("Phycicoccus")), expression(italic("Brevibacterium")), 
                        expression(italic("Agromyces")), expression(italic("Microbacterium")), 
                        expression(italic("Tetrasphaera")), expression(italic("Streptomyces"))),
         annotation_legend = TRUE)
dev.off()

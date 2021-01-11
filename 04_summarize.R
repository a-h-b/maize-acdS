# Analysis of acdS data
# Lucie Gebauer, last edited 12/2020


OTUs <- read.table("data/01_01_rarefied_otus_ed.txt")[,6:4623]
Treatments <- read.table("data/01_01_rarefied_otus_ed.txt")[,1:5]
phylum<-t(read.table("data/01_otu_ed.txt")[,62])
class<-t(read.table("data/01_otu_ed.txt")[,63])
order<-t(read.table("data/01_otu_ed.txt")[,64])
family<-t(read.table("data/01_otu_ed.txt")[,65])
genus<-t(read.table("data/01_otu_ed.txt")[,66])
species<-t(read.table("data/01_otu_ed.txt")[,67])

#check levels
levels(phylum)
levels(class)
levels(species)
levels(genus)
levels(order)
levels(family)

#---
Phylum_tab<-OTUs 
colnames(Phylum_tab)<-phylum
Phylum_table<-cbind(Treatments, Phylum_tab)

write.table(Phylum_table, "data/04_Phylum_table_rare.txt", sep="\t", col.names=TRUE, row.names=TRUE)

#---

Class_tab<-OTUs
colnames(Class_tab)<-class
Class_table<-cbind(Treatments, Class_tab)

write.table(Class_table, "data/04_Class_table_rare.txt", sep="\t", col.names=TRUE, row.names=TRUE)

#---

Order_tab<-OTUs
colnames(Order_tab)<-order
Order_table<-cbind(Treatments, Order_tab)

write.table(Order_table, "data/04_Order_table_rare.txt", sep="\t", col.names=TRUE, row.names=TRUE)

#---

Family_tab<-OTUs
colnames(Family_tab)<-family
Family_table<-cbind(Treatments, Family_tab)

write.table(Family_table, "data/04_Family_table_rare.txt", sep="\t", col.names=TRUE, row.names=TRUE)

#---

Genus_tab<-OTUs
colnames(Genus_tab)<-genus
Genus_table<-cbind(Treatments, Genus_tab)

write.table(Genus_table, "data/04_Genus_table_rare.txt", sep="\t", col.names=TRUE, row.names=TRUE)

#---

Species_tab<-OTUs
colnames(Species_tab)<-species
Species_table<-cbind(Treatments, Species_tab)

write.table(Species_table, "data/04_Species_table_rare.txt", sep="\t", col.names=TRUE, row.names=TRUE)

#---

precol <- Phylum_table[,1:6]

#1 Phylum_table_sum

names(Phylum_tab)<- gsub("\\..*", "", names(Phylum_tab))
unames<-unique(gsub("\\..*", "", names(Phylum_tab)))

P_summary<- sapply(unames, function(y) rowSums(Phylum_tab[y == names(Phylum_tab)]))
P_final <- cbind(Treatments, P_summary)

write.table(P_final, "data/04_Phylum_table_rare_sum.txt", sep = "\t", col.names=TRUE, row.names=TRUE)

#2 Class

names(Class_tab)<- gsub("\\..*", "", names(Class_tab))
unames<-unique(gsub("\\..*", "", names(Class_tab)))

C_summary<- sapply(unames, function(y) rowSums(Class_tab[y == names(Class_tab)]))
C_final <- cbind(Treatments, C_summary)

write.table(C_final, "data/04_Class_table_rare_sum.txt", sep = "\t", col.names=TRUE, row.names=TRUE)

#3 Order

names(Order_tab)<- gsub("\\..*", "", names(Order_tab))
unames<-unique(gsub("\\..*", "", names(Order_tab)))

O_summary<- sapply(unames, function(y) rowSums(Order_tab[y == names(Order_tab)]))
O_final <- cbind(Treatments, O_summary)

write.table(O_final, "data/04_Order_table_rare_sum.txt", sep = "\t", col.names=TRUE, row.names=TRUE)

#4 Family

names(Family_tab)<- gsub("\\..*", "", names(Family_tab))
unames<-unique(gsub("\\..*", "", names(Family_tab)))

F_summary<- sapply(unames, function(y) rowSums(Family_tab[y == names(Family_tab)]))
F_final <- cbind(Treatments, F_summary)

write.table(F_final, "data/04_Family_table_rare_sum.txt", sep = "\t", col.names=TRUE, row.names=TRUE)

#5 Genus

names(Genus_tab)<- gsub("\\..*", "", names(Genus_tab))
unames<-unique(gsub("\\..*", "", names(Genus_tab)))

G_summary<- sapply(unames, function(y) rowSums(Genus_tab[y == names(Genus_tab)]))
G_final <- cbind(Treatments, G_summary)

write.table(G_final, "data/04_Genus_table_rare_sum.txt", sep = "\t", col.names=TRUE, row.names=TRUE)

#6 Species

names(Species_tab)<- gsub("\\..*", "", names(Species_tab))
unames<-unique(gsub("\\..*", "", names(Species_tab)))

S_summary<- sapply(unames, function(y) rowSums(Species_tab[y == names(Species_tab)]))
S_final <- cbind(Treatments, S_summary)

write.table(S_final, "data/04_Species_table_rare_sum.txt", sep = "\t", col.names=TRUE, row.names=TRUE)


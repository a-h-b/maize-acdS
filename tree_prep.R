condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))
library(Biostrings)
library(vegan)
library(ape)

a <- read.delim("00_otu_table.corrected.tsv",stringsAsFactors=F)
f <- readDNAStringSet("01_01_rarefied_otus_ed-fromCorrected.fasta")
f2 <- readDNAStringSet("OTUs_Ref.mafft.fasttree.cluster98.fasta")
cl <- read.delim("OTUs_Ref.mafft.fasttree.cluster98.tsv")
tax <- read.delim("../../../DBs/amplicon/acdSgene.shortname.taxonomy",stringsAsFactors=F,header=F)
tr <- read.tree("OTUs_Ref.mafft.fasttree.cluster98.mafft.fasttree.newick")


tr2 <- drop.tip(tr,"OTU_009476")
write.tree(tr2,"OTUs_Ref.mafft.fasttree.cluster98.mafft.fasttree.drop1.newick")

dfr <- data.frame("name"=tr2$tip.label[!grepl("_",tr2$tip.label)],"tax"=sapply(tr2$tip.label[!grepl("_",tr2$tip.label)],function(x)tax$V2[tax$V1==x]),stringsAsFactors=F)
a2 <- a[ a$OTU %in% tr2$tip.label,]
dfr$leaf_dot_color <- "k_grey"
taxo <- gsub(" ","_",apply(a2[,c("kingdom","phylum","class","order","family","genus","species" )],1,function(x)paste(x,collapse=";",sep=";")))
a2$taxonomy <- taxo
dfo <- data.frame("name"=tr2$tip.label[grep("OTU",tr2$tip.label)],"tax"=sapply(tr2$tip.label[grep("OTU",tr2$tip.label)],function(x)a2$taxonomy[a2$OTU==x]),stringsAsFactors=F)
dfo$leaf_dot_color <- "k_yellow_green"
df <- rbind(dfo,dfr)
df$cluster <- sapply(df$name,function(x) cl$cluster[cl$OTU==x])
df$OTUs <- sapply(df$cluster,function(x) length(which(cl$cluster==x&grepl("OTU",cl$OTU))))
df$refs <- sapply(df$cluster,function(x) length(which(cl$cluster==x&!grepl("OTU",cl$OTU))))
df$bar1_height <- sqrt(df$refs)
df$bar2_height <- sqrt(df$OTUs)
df$bar1_color <- "k_grey"
df$bar2_color <- "k_yellow_green"
df$shortTax <- sapply(df$tax,function(x) paste(unlist(strsplit(x,split=";"))[1:4],collapse=";",sep=";"))
df$leaf_label_color <- ""
df$leaf_label_color[df$shortTax %in% c("unclassified;unclassified;unclassified;unclassified","Bacteria;;;")] <- "k_grey"
df$leaf_label_color[grep("Eukaryota",df$shortTax)] <- "k_buff"
df$leaf_label_color[df$shortTax %in% c("Bacteria;Actinobacteria;Actinobacteria;","Bacteria;Actinobacteria;Actinobacteria;unclassified_Actinobacteria")] <- "r_royalblue"
df$leaf_label_color[df$shortTax=="Bacteria;Actinobacteria;Actinobacteria;Corynebacteriales" ] <- "r_turquoise1"
df$leaf_label_color[df$shortTax=="Bacteria;Actinobacteria;Actinobacteria;Geodermatophilales" ] <- "r_turquoise3"
df$leaf_label_color[df$shortTax=="Bacteria;Actinobacteria;Actinobacteria;Micrococcales" ] <- "r_turquoise4"
df$leaf_label_color[df$shortTax=="Bacteria;Actinobacteria;Actinobacteria;Micromonosporales" ] <- "r_skyblu2"
df$leaf_label_color[df$shortTax=="Bacteria;Actinobacteria;Actinobacteria;Micromonosporales" ] <- "r_skyblue2"
df$leaf_label_color[df$shortTax=="Bacteria;Actinobacteria;Actinobacteria;Nakamurellales" ] <- "r_skyblue4"
df$leaf_label_color[df$shortTax=="Bacteria;Actinobacteria;Actinobacteria;Propionibacteriales" ] <- "r_slateblue2"
df$leaf_label_color[df$shortTax=="Bacteria;Actinobacteria;Actinobacteria;Pseudonocardiales" ] <- "r_slateblue4"
df$leaf_label_color[df$shortTax=="Bacteria;Actinobacteria;Actinobacteria;Streptomycetales" ] <- "r_navy"
df$leaf_label_color[df$shortTax=="Bacteria;Proteobacteria;;" ] <- "r_goldenrod"
df$leaf_label_color[df$shortTax=="Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales" ] <- "r_sandybrown"
df$leaf_label_color[df$shortTax=="Bacteria;Proteobacteria;Alphaproteobacteria;Rhodobacterales" ] <- "r_sienna2"
df$leaf_label_color[df$shortTax=="Bacteria;Proteobacteria;Alphaproteobacteria;Rhodospirillales" ] <- "r_sienna4"
df$leaf_label_color[df$shortTax=="Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales" ] <- "r_yellow4"
df$leaf_label_color[df$shortTax=="Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales" ] <- "r_olivedrab2"
df$leaf_label_color[df$shortTax=="Bacteria;Proteobacteria;Gammaproteobacteria;Oceanospirillales" ] <- "r_limegreen"
df$leaf_label_color[df$shortTax=="Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales" ] <- "r_lawngreen"
df$leaf_label_color[df$shortTax=="Bacteria;Proteobacteria;Gammaproteobacteria;unclassified_Gammaproteobacteria" ] <- "r_darkolivegreen"
df$new_name <- paste(df$name,"_",gsub(".+;","",gsub(";+$","",df$tax)),sep="")

write.table(df[,c("name","leaf_dot_color","bar1_height","bar2_height","bar1_color","bar2_color","leaf_label_color","new_name")],"mapping.txt",sep="\t",quote=F,row.names=F)

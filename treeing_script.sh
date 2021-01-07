cat 01_01_rarefied_otus_ed-fromCorrected.fasta ../References.fasta >> OTUs_Ref.fasta
mafft OTUs_Ref.fasta > OTUs_Ref.mafft.fasta
FastTreeMP -nt -gamma -no2nd -fastest -spr 4 -log fasttree.log -quiet OTUs_Ref.mafft.fasta > OTUs_Ref.mafft.fasttree.newick


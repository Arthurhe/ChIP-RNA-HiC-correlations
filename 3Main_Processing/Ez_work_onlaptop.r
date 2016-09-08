#pick target
target="K562"
interactionlist_pos=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_center_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_negative_loops.txt",stringsAsFactors = F,sep="\t")

target="GM12"
interactionlist_pos=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/GM12878_25K_center_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("",stringsAsFactors = F,sep="\t")

target_markers=c("RNA","CAGE","DNase","FAIRE","RRBS","ATF3","BCL3","BCLAF1","BHLHE40","CEBPB","CHD1","CHD2","CREB1","CTCF","CUX1",
                 "E2F4","EGR1","ELF1","ELK1","EP300","ETS1","EZH2","FOS","GABPA","H2AZ","H3K27ac","H3K27me3","H3K36me3","H3K4me1",
                 "H3K4me2","H3K4me3","H3K79me2","H3K9ac","H3K9me3","H4K20me1","JUND","MAFK","MAX","MAZ","MEF2A","MYC","NFE2","NFYA",
                 "NFYB","NR2C2","NRF1","PML","POLR2A","POLR3G","RAD21","RCOR1","REST","RFX5","SIX5","SMC3","SPI1","SP1","SRF","STAT1",
                 "STAT5A","TAF1","TBL1XR1","TBP","USF1","USF2","YY1","ZBTB33","ZNF143","ZNF274","ZNF384","HCFC1")

#loop length check (sampler quality check)
par(mfrow=c(2,1))
plot(density(interactionlist_pos[,4],bw = 10000),ylim=c(0,5e-7),main="loop length distribution of real interactions",
     ylab="probabilty",xlab="loop length",sub=paste("totally",nrow(interactionlist_pos),"samples"))
plot(density(interactionlist_neg[,4],bw = 10000),ylim=c(0,5e-7),main="loop length distribution of fake interactions",
     ylab="probabilty",xlab="loop length",sub=paste("totally",nrow(interactionlist_neg),"samples"))

#region interactivities statistics (all location have interactions)
active_sites=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_12.5Kres_active_unique_anchors.bed",stringsAsFactors = F,sep="\t")
negative_sites=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_12.5Kres_negative_anchors.bed",stringsAsFactors = F,sep="\t")

plot(hist(c(active_sites[,4],negative_sites[,4]),seq(-0.5,199.5,1),freq=F),xlim=c(0,50),type="l",lwd=2,ylim=c(0,0.12),
     main="25K regions interactivities",ylab="percentage",xlab="num of region interact with")

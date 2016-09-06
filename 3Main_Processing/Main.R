#Interaction statistics
read the shit
length distribution
bin not overlap with target


#Sampler
interactionlist=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/GM12878_25K_center_interactions_simplified.txt",stringsAsFactors = F,sep="\t")

interactionlist_pos=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_center_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_negative_loops.txt",stringsAsFactors = F,sep="\t")

par(mfrow=c(2,1))
plot(density(interactionlist_pos[,4],bw = 10000),ylim=c(0,5e-7),main="loop length distribution of real interactions",
     ylab="probabilty",xlab="loop length",sub=paste("totally",nrow(interactionlist_pos),"samples"))
plot(density(interactionlist_neg[,4],bw = 10000),ylim=c(0,5e-7),main="loop length distribution of fake interactions",
     ylab="probabilty",xlab="loop length",sub=paste("totally",nrow(interactionlist_neg),"samples"))

#all location have interactions
active_sites=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_12.5Kres_active_unique_anchors.bed",stringsAsFactors = F,sep="\t")
negative_sites=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_12.5Kres_negative_anchors.bed",stringsAsFactors = F,sep="\t")

plot(hist(c(active_sites[,4],negative_sites[,4]),seq(-0.5,199.5,1),freq=F),xlim=c(0,50),type="l",lwd=2,ylim=c(0,0.12),
          main="25K regions interactivities",ylab="percentage",xlab="num of region interact with")
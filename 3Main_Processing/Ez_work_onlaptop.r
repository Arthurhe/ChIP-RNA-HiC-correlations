library(gplots)
library(dtw)
library(Rtsne) ##Rtsne 
source("F:/DATA/R/CTCFanalysis_functionPack.R")

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
pdf(file="loop length distribution.pdf",width=5,height=8)
par(mfrow=c(2,1))
plot(density(interactionlist_pos[,4],bw = 10000),ylim=c(0,5e-7),main="loop length distribution of real interactions",
     ylab="probabilty",xlab="loop length",sub=paste("totally",nrow(interactionlist_pos),"samples"))
plot(density(interactionlist_neg[,4],bw = 10000),ylim=c(0,5e-7),main="loop length distribution of fake interactions",
     ylab="probabilty",xlab="loop length",sub=paste("totally",nrow(interactionlist_neg),"samples"))
dev.off()

#region interactivities statistics (all location have interactions)
active_sites=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_12.5Kres_active_unique_anchors.bed",stringsAsFactors = F,sep="\t")
negative_sites=read.table("F:/DATA/R/Kees/1608_HicChipRNACor/data/HiC/K562_25K_12.5Kres_negative_anchors.bed",stringsAsFactors = F,sep="\t")

plot(hist(c(active_sites[,4],negative_sites[,4]),seq(-0.5,199.5,1),freq=F),xlim=c(0,50),type="l",lwd=2,ylim=c(0,0.12),
     main="25K regions interactivities",ylab="percentage",xlab="num of region interact with")

matrixlist_500=c()
newmatrixlist=c()
characteristic_marker=read.table("file:///F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/characteristic_marker_DMD.txt",sep="\t")
for(i in 1:50){
  matrixlist_500[[i]]=read.csv(paste("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/genome_looping_matrix_",i,sep=""))[,-1]  
  matrixlist_500[[i]][,"RRBS"]=matrixlist_500[[i]][,"RRBS"]/100
  if(sum(characteristic_marker[,i])!=0){newmatrixlist[[i]]=rowMeans(matrixlist_500[[i]][,characteristic_marker[,i]])}
  else{newmatrixlist[[i]]=rowMeans(matrixlist_500[[i]])}
}

#check DMD
importance=matrix(0,71,2)
rownames(importance)=target_markers
for(i in 1:50){
  temp=unique(characteristic_marker[,i])
  for(i in temp){
    importance[i,2]=importance[i,2]+1}
}

#target_markers[order(importance[,2],decreasing=T)]
#importance[order(importance[,2],decreasing=T),2]
par(mfrow=c(1,3))
linecol=c(rep("red",25),rep("black",25))
plot(0,type="n",ylim=c(0,25.5),xlim=c(0,max(sapply(1:length(matrixlist_500),function(x){nrow(matrixlist_500[[x]])}))),main="average signal\noriginal")
for(i in 1:length(matrixlist_500)){
  lines(cbind(1:nrow(matrixlist_500[[i]]),rowMeans(matrixlist_500[[i]])+i*0.5-1),col=linecol[i])
}
plot(0,type="n",ylim=c(0,25.5),xlim=c(0,max(sapply(1:length(matrixlist_500),function(x){nrow(matrixlist_500[[x]])}))),main="average signal\nchanged by DMD")
for(i in 1:length(matrixlist_500)){
  lines(cbind(1:nrow(matrixlist_500[[i]]),(newmatrixlist[[i]]-rowMeans(matrixlist_500[[i]]))+i*0.5-1),col=linecol[i])
}
plot(0,type="n",ylim=c(0,25.5),xlim=c(0,max(sapply(1:length(matrixlist_500),function(x){nrow(matrixlist_500[[x]])}))),main="average signal\nDMD processed")
for(i in 1:length(matrixlist_500)){
  lines(cbind(1:nrow(matrixlist_500[[i]]),newmatrixlist[[i]]+i*0.5-1),col=linecol[i])
}

par(mfrow=c(1,2))
dis_matrix_500=dtw_distMatCreat(matrixlist_500)
rtsne_out2=Rtsne(dis_matrix_500,dims=2,is_distance=T,perplexity=10)
plot(rtsne_out2$Y,col=linecol,pch=16,main="DTW_original")
dis_matrix_500=dtw_distMatCreat(newmatrixlist)
rtsne_out2=Rtsne(dis_matrix_500,dims=2,is_distance=T,perplexity=10)
plot(rtsne_out2$Y,col=linecol,pch=16,main="DTW_DMD_processed")

#testing on 200 samples.
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/matrix_200_short.rdata")
load("~/Analysis/201608_HicChipRnaCor/data/matrix_200_short.rdata")
learning_list_200=matrixlist_500

#initiate all necessary variable
featurelib=c()
learning_list_200_feature_vec=c()
row_vec=c()
time_vec=c()

#create feature list
for(i in 1:length(learning_list_200)){
  ptm <- proc.time()
  learning_list_200_feature_vec[[i]]=MatrixScan_featureLibBuild(learning_list_200[[i]],wd=1,distanceThreshold=0.2)
  row_vec[i]=nrow(learning_list_200[[i]])
  time_vec[i]=(proc.time()-ptm)[3]
  print(paste(i,"|",row_vec[i],"rows |",time_vec[i]))
}

#examing the result of 200 samples:
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/featureAnaOn200Short.Rdata")

featurelib_matrix=do.call(rbind,featurelib)
d=dist(featurelib_matrix,method = "euclidean")
fit=hclust(d) 
cluster_idx <- cutree(fit, 100)

featureMatrix=FeatureMatrixGenerator(learning_list_200_feature_vec,featurelib,cluster_idx)
idf_vec=IDF_calculator(featureMatrix)
tf_matrix=TF_calculator(featureMatrix)
tf_idf_matrix=TF_IDF_calculator(tf_matrix,idf_vec)
tf_idf_matrix=tf_idf_matrix[rowSums(tf_idf_matrix)!=0,colSums(tf_idf_matrix)!=0]


pca=prcomp(tf_idf_matrix)
plot(pca$x[,1:2],col=c(rep("red",100),rep("black",100)))

rtsne_out2=Rtsne(tf_idf_matrix,dims=2,perplexity=10)
plot(rtsne_out2$Y,col=c(rep("red",100),rep("black",100)))

heatmap.2(tf_idf_matrix,breaks=seq(0,0.5,0.025),
          col=colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")),
          trace="none",density.info="histogram",Colv=F,Rowv=F,notecol="black",dendrogram = "none",
          lhei=c(1,3),lwid=c(1,3),margins = c(5,5))

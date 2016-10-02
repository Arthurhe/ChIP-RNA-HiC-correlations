library(gplots)
library(ggplot2)
#library(dtw)
#library(Rtsne) ##Rtsne 
library(matrixStats)
library(e1071)
library(proxy) #dist
library(fastcluster)
library(cluster)
source("F:/DATA/R/Kees/1608_HicChipRNACor/3Main_Processing/HicChipRNACor_functionPack.r")

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

#################################################
#examing the result of 15000 samples:
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/featureOn15000Short.Rdata")
learning_feature_vec[2001:10000]=scanning_feature_vec[2001:10000]
#find unique feature on looping region #filter 1
alltruefeature=do.call(c,learning_feature_vec[1001:5000])
allfalsefeature=do.call(c,learning_feature_vec[5001:9000])
feature_stat_pos=Feature_numstat_vec(alltruefeature,nrow(featurelib_matrix))
feature_stat_neg=Feature_numstat_vec(allfalsefeature,nrow(featurelib_matrix))
uniqueFeature=which(abs(log10((feature_stat_pos)/(feature_stat_neg)))>=1 & (feature_stat_pos+feature_stat_neg)>=5)

#adding none unique feature
#calculate the important centers
d=dist(featurelib_matrix[uniqueFeature,],method = "manhattan",upper = T) 
fit=hclust(d,method="complete") 
rm(d)
cluster_idx <- as.vector(cutree(fit, h=0.3))
#calculate the cluster_centers
uniqueFeatureCenters=GetGroupCenter(featurelib_matrix[uniqueFeature,],cluster_idx)
#feature assign to centers
d=dist(featurelib_matrix,uniqueFeatureCenters,method = "manhattan")
cluster_idx=apply(d,1,which.min)*(rowMins(d)<=0.25)
cluster_idx[cluster_idx==0]=max(cluster_idx)+1
#sum(cluster_idx!=max(cluster_idx))

#only uniqueFeature #best
cluster_idx=rep(0,nrow(featurelib_matrix))
cluster_idx[uniqueFeature]=1:length(uniqueFeature)
cluster_idx[cluster_idx==0]=max(cluster_idx)+1

#positional calculation
posfeatureDistriList=FeatureDistriListGenrator_abs(learning_feature_vec[1001:5000],cluster_idx)
negfeatureDistriList=FeatureDistriListGenrator_abs(learning_feature_vec[5001:9000],cluster_idx)

featureMatrix=FeaturePositionScoreMatrixGenerator(c(learning_feature_vec,K562_scanning_feature_vec),cluster_idx,featureDistriList)
featureMatrix=featureMatrix[,-max(cluster_idx)]

#feature score calculation for each sample
#clustering=clara(featurelib_matrix[trueUniquefaeture,], 1000, metric = "manhattan", stand=F, sample=50, medoids.x = F)
featureMatrix=FeatureMatrixGenerator(c(learning_feature_vec,K562_scanning_feature_vec),nrow(featurelib_matrix),cluster_idx) #clustering$clustering
cluster_specificity=log2(colSums(featureMatrix[1001:5000,],na.rm = T)/colSums(featureMatrix[5001:9000,],na.rm = T))
ignore=which(abs(cluster_specificity)<2)
tf_matrix=TF_calculator(featureMatrix)

#idf is useless here
idf_vec=IDF_calculator(featureMatrix)
tf_idf_matrix=TF_IDF_calculator(tf_matrix,idf_vec)
featureMatrix=tf_idf_matrix

#svm start
temp=svm(featureMatrix[1001:9000,],c(rep("1",4000),rep("0",4000)),type="C-classification",cost=1000,scale=F,probability = T)
prediction <- predict(temp, featureMatrix, probability = T)
probability=attr(prediction, "probabilities")[,1]
plot(probability,pch=19,cex=0.6)
prediction = as.numeric(levels(prediction))[prediction]

nas=which(is.na(rowSums(featureMatrix)))
breakpoint=c(0,1000,5000,9000,10000,12500,15000)
for(i in 1:(length(breakpoint)-1)){
  predictedTure=sum(prediction[(breakpoint[i]+1-sum(nas<breakpoint[i])):(breakpoint[i+1]-sum(nas<breakpoint[i+1]))])
  print(paste(predictedTure,breakpoint[i+1]-sum(nas<breakpoint[i+1])-breakpoint[i]+sum(nas<breakpoint[i])-predictedTure,breakpoint[i+1]-sum(nas<breakpoint[i+1])-breakpoint[i]+sum(nas<breakpoint[i])))
}

plot(sort(log2(colMeans(featureMatrix[1:4000,],na.rm = T)-colMeans(featureMatrix[4001:8000,],na.rm = T))))






temp=svm(tf_idf_matrix[1:8000,],c(rep(1,4000),rep(0,4000)),cost=10000,scale=F)
prediction <- predict(temp, tf_idf_matrix)
plot(prediction,pch=19)

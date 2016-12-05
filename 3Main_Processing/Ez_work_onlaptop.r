library(gplots)
library(ggplot2)
#library(dtw)
#library(Rtsne) ##Rtsne 
library(matrixStats)
library(e1071)
library(proxy) #dist
library(fastcluster)
library(cluster)
library(ROCR)
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
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/featureOnAllGM_short.Rdata")
rm("learning_feature_vec")

#find unique feature on looping region #filter 1
alltruefeature=do.call(c,scanning_feature_vec[27399:32398])
allfalsefeature=do.call(c,scanning_feature_vec[32399:37398])
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
#feature assign to centers(uniqueFeature) #merge 500 more feature, given 2968 unique feature
d=dist(featurelib_matrix,featurelib_matrix[uniqueFeature,],method = "manhattan")
cluster_idx=apply(d,1,which.min)*(rowMins(d)<=0.25)
cluster_idx[cluster_idx==0]=max(cluster_idx)+1
rm(d)
peak_region_matrix=NULL

#feature normalizing to all 1s
featurelib_matrix[featurelib_matrix>0.05]=1
featurelib_matrix[featurelib_matrix<=0.05]=0
uniqueImportantFeature=featurelib_matrix[uniqueFeature[!duplicated(featurelib_matrix[uniqueFeature,])],]
d=dist(featurelib_matrix,uniqueImportantFeature,method = "manhattan")
cluster_idx=apply(d,1,which.min)*(rowMins(d)<=1)
cluster_idx[cluster_idx==0]=max(cluster_idx)+1
peak_region_matrix=NULL

#positional calculation
posfeatureDistriList=FeatureDistriListGenrator_abs(learning_feature_vec[1001:5000],cluster_idx)
negfeatureDistriList=FeatureDistriListGenrator_abs(learning_feature_vec[5001:9000],cluster_idx)
peak_region_matrix=FeaturePositionScoreMatrixGenerator(posfeatureDistriList,negfeatureDistriList,0.001)

#feature score calculation for each sample
#clustering=clara(featurelib_matrix[trueUniquefaeture,], 1000, metric = "manhattan", stand=F, sample=50, medoids.x = F)
featureMatrix=FeatureMatrixGenerator(c(learning_feature_vec,K562_scanning_feature_vec),nrow(featurelib_matrix),cluster_idx,peak_region_matrix) #clustering$clustering
cluster_specificity=log2(colSums(featureMatrix[1001:5000,],na.rm = T)/colSums(featureMatrix[5001:9000,],na.rm = T))
#ignore=which(abs(cluster_specificity)<2)
tf_matrix=TF_calculator(featureMatrix)
featureMatrix=tf_matrix

#idf is useless here
idf_vec=IDF_calculator(featureMatrix)
tf_idf_matrix=TF_IDF_calculator(tf_matrix,idf_vec)
featureMatrix=tf_idf_matrix

#svm start
nas=which(is.na(rowSums(featureMatrix)))
featureMatrix[nas,]=0
svm_model=svm(featureMatrix[2001:8000,],c(rep("1",3000),rep("0",3000)),type="C-classification",cost=1000,scale=F,probability = T)
pred <- predict(svm_model, featureMatrix, probability = T)
prob=attr(pred, "probabilities")[,1]
pred=as.numeric(levels(pred))[pred]
#re-arrange pred and prob
fill_seq=c(3,1,2,4,5,6)
breakpoint=c(0,2000,5000,8000,10000,12500,15000)
lb=as.numeric(c(rep(1,5000),rep(0,5000),rep(1,2500),rep(0,2500)))
prob_arranged=c()
pred_arranged=c()
lb_arranged=c()
breakpoint_arranged=rep(0,length(breakpoint))
for(i in 1:(length(breakpoint)-1)){
  start=breakpoint[(which(fill_seq==i))]+1
  end=breakpoint[(which(fill_seq==i))+1]
  prob_arranged=c(prob_arranged,prob[start:end])
  pred_arranged=c(pred_arranged,pred[start:end])
  breakpoint_arranged[i+1]=breakpoint_arranged[i]+end-start+1
  lb_arranged=c(lb_arranged,lb[start:end])
}

#calculate performance
performance_table=matrix(0,7,3)
for(i in 1:(length(breakpoint)-1)){
  predictedTure=sum(pred_arranged[(breakpoint_arranged[i]+1):(breakpoint_arranged[i+1])],na.rm = T)
  performance_table[i,1]=predictedTure
  performance_table[i,2]=breakpoint_arranged[i+1]-breakpoint_arranged[i]-predictedTure
  performance_table[i,3]=breakpoint_arranged[i+1]-breakpoint_arranged[i]
}
performance_table[7,]=colSums(performance_table[1:6,])
precision=round(c(performance_table[1,1]/sum(performance_table[1:2,1]),performance_table[3,1]/sum(performance_table[3:4,1]),performance_table[5,1]/sum(performance_table[5:6,1])),2)
sensitivity=round(c(performance_table[1,1]/sum(performance_table[1,1:2]),performance_table[3,1]/sum(performance_table[3,1:2]),performance_table[5,1]/sum(performance_table[5,1:2])),2)

par(mfrow=c(1,3))
plot(prob_arranged[1:6000],pch=19,cex=0.6,main="training set",ylab="predicted looping probability",xlab="sample")
abline(v=3000.5,col="red",lwd=2)
plot(prob_arranged[6001:10000],pch=19,cex=0.6,main="test set - same celltype",ylab="predicted looping probability",xlab="sample")
abline(v=2000.5,col="red",lwd=2)
plot(prob_arranged[10001:15000],pch=19,cex=0.6,main="test set - other celltype",ylab="predicted looping probability",xlab="sample")
abline(v=2500.5,col="red",lwd=2)


main_title_list=c("training set","test set - same celltype","test set - other celltype")
collist=c("dodgerblue3","gold2","firebrick3")
ROC_PRC=matrix(0,2,3)
plotlist=c("AUC of ROC","AUC of PRC")
par(mfrow=c(2,1))
for(k in 1:2){
  for(i in 1:((length(breakpoint)-1)/2)){
    ROCR_pred=ROCR::prediction(prob_arranged[(breakpoint_arranged[i*2-1]+1):breakpoint_arranged[1+i*2]],
                               lb_arranged[(breakpoint_arranged[i*2-1]+1):breakpoint_arranged[1+i*2]])
    if(k==1){
      perf=ROCR::performance(ROCR_pred,"tpr","fpr")
      ROC_PRC[k,i]=PRROC::roc.curve(scores.class0 = prob_arranged[(breakpoint_arranged[i*2-1]+1):breakpoint_arranged[i*2]],scores.class1 = prob_arranged[(breakpoint_arranged[i*2]+1):breakpoint_arranged[1+i*2]])$auc
    }else{
      perf=ROCR::performance(ROCR_pred,"prec","rec")
      ROC_PRC[k,i]=PRROC::pr.curve(scores.class0 = prob_arranged[(breakpoint_arranged[i*2-1]+1):breakpoint_arranged[i*2]],scores.class1 = prob_arranged[(breakpoint_arranged[i*2]+1):breakpoint_arranged[1+i*2]])$auc.integral
    }
    if(i==1){plot(perf,col=collist[i],lwd=2,ylim=c(0,1))
    }else{plot(perf,col=collist[i],add = TRUE,lwd=2)}
  }
  title(plotlist[k])
  legend("bottomright",paste(main_title_list," (",round(ROC_PRC[k,],2),")",sep=""),bty = "n",lty=c(1,1),lwd=c(2,2),col=collist)
}





plot(perf,colorize=TRUE)
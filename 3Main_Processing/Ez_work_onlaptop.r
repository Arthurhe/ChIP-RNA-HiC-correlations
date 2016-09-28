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
#feature build_scan_analysis
#build feature set in sampleset 200_short_1
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/matrix_200_short.rdata")
load("~/Analysis/201608_HicChipRnaCor/data/matrix_200_short.rdata")
learning_list_200=matrixlist_500

#initiate all necessary variable
learning_list_200_feature_vec=c()
row_vec=c()
time_vec=c()
featurelib_matrix=NA
cluster_info=list("cluster"=0)
#create feature list
for(i in 1:length(learning_list_200)){
  ptm <- proc.time()
  list[learning_list_200_feature_vec[[i]],featurelib_matrix,cluster_info]=MatrixScan_featureLibBuild_Advance(learning_list_200[[i]],wd=1,featurelib_matrix,cluster_info,distanceThreshold=0.2)
  row_vec[i]=nrow(learning_list_200[[i]])
  time_vec[i]=(proc.time()-ptm)[3]
  print(paste(i,"|",row_vec[i],"rows |",time_vec[i]))
}
save(featurelib,learning_list_200_feature_vec,row_vec,time_vec,file="featureAnaOn200Short.Rdata")

#scan feature set on sampleset 200_short_2
load("~/Analysis/201608_HicChipRnaCor/data/matrix_200_short_2.rdata")
load("~/Analysis/201608_HicChipRnaCor/data/featureAnaOn200Short.Rdata")
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/featureAnaOn200Short.Rdata")
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/matrix_200_short_2.rdata")

#initiate all necessary variable
scanning_list_200_feature_vec=c()
scan_row_vec=c()
scan_time_vec=c()
#featurelib_matrix=do.call(rbind,featurelib)
#cluster_info=KmeanFeatureClustering(featurelib_matrix,as.integer(nrow(featurelib_matrix)/250))
for(i in 1:length(matrixlist_500)){
  ptm <- proc.time()
  scanning_list_200_feature_vec[[i]]=MatrixScan_Advance2(matrixlist_500[[i]],featurelib_matrix,cluster_info,distanceThreshold=0.2)
  scan_row_vec[i]=nrow(matrixlist_500[[i]])
  scan_time_vec[i]=(proc.time()-ptm)[3]
  print(paste(i,"|",scan_row_vec[i],"rows |",scan_time_vec[i]))
}
save(featurelib,scanning_list_200_feature_vec,scan_row_vec,scan_time_vec,file="featureScanOn200Short2.Rdata")

for(i in 1:length(scanning_list_200_feature_vec)){if(sum(scanning_list_200_feature_vec[[i]]!=learning_list_200_feature_vec[[i]])!=0){print(i)}}


#examing the result of 200+200 samples:
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/featureAnaOn200Short.Rdata")
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/featureScanOn200Short2.Rdata")

#featurelib_matrix=do.call(rbind,featurelib)
#featurelib_matrix=apply(featurelib_matrix,2,function(x){x/mean(x)})
d=dist(featurelib_matrix,method = "manhattan",upper = T) #euclidean
#system.time(d<-dist(featurelib_matrix,method = "manhattan",upper = T))
fit=hclust(d,method="complete") 
cluster_idx <- as.vector(cutree(fit, 100))
#c(learning_list_200_feature_vec,scanning_list_200_feature_vec)
featureMatrix=FeatureMatrixGenerator(learning_list_200_feature_vec,nrow(featurelib_matrix),cluster_idx)
idf_vec=IDF_calculator(featureMatrix)
tf_matrix=TF_calculator(featureMatrix)
tf_idf_matrix=TF_IDF_calculator(tf_matrix,idf_vec)
#tf_idf_matrix=tf_idf_matrix[rowSums(tf_idf_matrix)!=0,colSums(tf_idf_matrix)!=0]

#visualize the feature belonging
heatmap.2(tf_idf_matrix,breaks=seq(0,0.5,0.025),
          col=colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")),
          trace="none",density.info="histogram",Colv=F,Rowv=F,notecol="black",dendrogram = "none",
          lhei=c(1,3),lwid=c(1,3),margins = c(5,5))

temp=cbind(tf_idf_matrix[,colSums(tf_idf_matrix[1:100,])-colSums(tf_idf_matrix[101:200,])>0.1],
           tf_idf_matrix[,colSums(tf_idf_matrix[1:100,])-colSums(tf_idf_matrix[101:200,])<(-0.1)])

heatmap.2(as.matrix(temp[1:200,]),
          col=colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")),
          trace="none",density.info="histogram",Colv=F,Rowv=F,notecol="black",dendrogram = "none",
          lhei=c(1,3),lwid=c(1,3),margins = c(5,5))

#visualize feature
k=c()
for(i in 1:max(cluster_idx)){k[i]=sum(cluster_idx==i)}
p <- ggplot(data = k, aes(x = factor(gene,levels=c("H3f3a","H3f3b","H3f3c","H3f3a_optimized")), y = mean,
                          fill = factor(celltype,levels=c("WT","KO","OE"))))
p + geom_bar(stat = "identity",position = position_dodge(0.9)) +  scale_fill_brewer(palette="Set1",name = "celltype")+
  labs(x = "gene", y = "normalized RNA-seq count") +
  ggtitle("Gene expression level of H3f3 genes")

temp=featurelib_matrix[cluster_idx==5,]
heatmap.2(as.matrix(temp),
          col=colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")),
          trace="none",density.info="histogram",Colv=F,Rowv=F,notecol="black",dendrogram = "none",
          lhei=c(1,3),lwid=c(1,3),margins = c(5,5))



temp=svm(tf_idf_matrix,c(rep(1,500),rep(0,500)),type="C-classification",cost=1000,scale=T)
prediction <- predict(temp, tf_idf_matrix)
prediction = as.numeric(levels(prediction))[prediction]
plot(prediction)
sum(prediction[c(201:300)]>0.5)
sum(prediction[c(301:400)]<0.5,na.rm = T)
tuned=tune(svm,tf_idf_matrix[1:200,],c(rep(1,100),rep(0,100)),ranges = list(cost=c(0.001,0.01,0.1,0,1,10,100,1000)))

#test on kmean

temp=kmeans(featurelib_matrix,10,iter.max = 100)
distance=dist(featurelib_matrix,temp$centers,method = "euclidean")
max(distance[(temp$cluster==1),1])
min(distance[(temp$cluster!=1),1])


#examing the result of 1000+1000+1000 samples:
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/uploading/featureOn11000Short.Rdata")

#find unique feature on looping region
alltruefeature=do.call(c,extrascanning_feature_vec[1:4000])
alltruefeature=unique(alltruefeature)
allfalsefeature=do.call(c,extrascanning_feature_vec[4001:8000])
allfalsefeature=unique(allfalsefeature)
trueUniquefaeture=union(setdiff(alltruefeature,allfalsefeature),setdiff(allfalsefeature,alltruefeature))

d=dist(featurelib_matrix[trueUniquefaeture,],method = "manhattan",upper = T) 
fit=hclust(d,method="complete") 
cluster_idx <- as.vector(cutree(fit, h=1))
clustering=rep(max(cluster_idx)+1,nrow(featurelib_matrix))
clustering[trueUniquefaeture]=cluster_idx
#clustering=clara(featurelib_matrix[trueUniquefaeture,], 1000, metric = "manhattan", stand=F, sample=50, medoids.x = F)
featureMatrix=FeatureMatrixGenerator(c(extrascanning_feature_vec,learning_feature_vec,scanning_feature_vec,K562_scanning_feature_vec),nrow(featurelib_matrix),clustering) #clustering$clustering
idf_vec=IDF_calculator(featureMatrix)
tf_matrix=TF_calculator(featureMatrix)
tf_idf_matrix=TF_IDF_calculator(tf_matrix,idf_vec)

temp=svm(tf_idf_matrix[1:8000,],c(rep(1,4000),rep(0,4000)),type="C-classification",cost=10000,scale=F)
prediction <- predict(temp, tf_idf_matrix)
prediction = as.numeric(levels(prediction))[prediction]
#plot(prediction,pch=19)

nas=which(is.na(rowSums(tf_idf_matrix)))
breakpoint=c(0,4000,8000,8500,9000,9500,10000,10500,11000)
for(i in 1:(length(breakpoint)-1)){
  predictedTure=sum(prediction[(breakpoint[i]+1-sum(nas<breakpoint[i])):(breakpoint[i+1]-sum(nas<breakpoint[i+1]))])
  print(paste(predictedTure,breakpoint[i+1]-sum(nas<breakpoint[i+1])-breakpoint[i]+sum(nas<breakpoint[i])-predictedTure,breakpoint[i+1]-sum(nas<breakpoint[i+1])-breakpoint[i]+sum(nas<breakpoint[i])))
}

MatrixScan_Advance1=function(regionMatrix,featurelib_matrix,cluster_info,distanceThreshold=0.2){
  #MatrixScan_cluster and base on similarity of last feature
  #regionMatrix: each row is a position on genome, each col is a marker
  #featurelib_matrix is a matrix feature that showed up before
  #cluster_info is the kmean object returned by kmean featurelib_matrix
  wd=ncol(featurelib_matrix)/71 #wd: the length of feature window
  if(wd%%1!=0){stop('incorrect featurelib_matrix col num')}
  if(!(length(featurelib_matrix)>1 & length(cluster_info)>1)){stop('featurelib_matrix & cluster_info must exist')}
  wdm=wd-1
  featureInRegion=rep(0,(nrow(regionMatrix)-wdm))
  #the first region
  feature=as.vector(regionMatrix[1:(1+wdm),])
  featureInRegion[1]=FindFeature_Advance1(feature,featurelib_matrix,cluster_info,distanceThreshold)
  #the following region
  for(x in 2:(nrow(regionMatrix)-wdm)){
    lastfeature=as.vector(regionMatrix[(x-1):(x+wdm-1),])
    feature=regionMatrix[x:(x+wdm),]
    if(sum(abs(lastfeature-feature))<distanceThreshold){featureInRegion[x]=featureInRegion[x-1]}
    else{featureInRegion[x]=FindFeature_Advance1(feature,featurelib_matrix,cluster_info,distanceThreshold)}
  }
  return(featureInRegion)
}

FindFeature_Advance1=function(feature,featurelib_matrix=NA,cluster_info=NA,distanceThreshold){
  #feature is a vector of one feature
  if(length(cluster_info)>1){
    distance=dist(feature,cluster_info$centers,method = "manhattan")
    cluster_idx=which.min(distance)
    pick_feature=which(cluster_info$cluster==cluster_idx)
    distance=dist(feature,featurelib_matrix[pick_feature,],method = "manhattan")
    if(min(distance)<distanceThreshold){tagFeature=pick_feature[which.min(distance)];return(tagFeature)}
    else{return(0)}
  }
  else if(length(featurelib)>1){
    distance=dist(feature,featurelib_matrix,method = "manhattan")
    if(distance[tagFeature]<distanceThreshold){tagFeature=which.min(distance);return(tagFeature)}
    else{return(0)}
  }
  else{
    return(0)
  }
}


MatrixScan_featureLibBuild=function(regionMatrix,wd=3,distanceThreshold=0.2){
  #regionMatrix: each row is a position on genome, each col is a marker
  #wd: the length of feature window
  #featurelib is a list of feature that showed up before
  #faeturelib must in the parent enviorment
  wdm=wd-1
  #the first region
  feature=regionMatrix[1:(1+wdm),]
  featureInRegion=rep(0,(nrow(regionMatrix)-wdm))
  featureInRegion[1]=FindFeature(feature,featurelib,distanceThreshold)
  if(featureInRegion[1]==0){
    featurelib_length=length(featurelib)+1
    featureInRegion[1]=featurelib_length
    featurelib[[featurelib_length]]<<-feature
  }
  #the following region
  for(x in 2:(nrow(regionMatrix)-wdm)){
    lastfeature=regionMatrix[(x-1):(x+wdm-1),]
    feature=regionMatrix[x:(x+wdm),]
    if(sum(abs(lastfeature-feature))<distanceThreshold){featureInRegion[x]=featureInRegion[x-1]}
    else{
      featureInRegion[x]=FindFeature(feature,featurelib,distanceThreshold)
      if(featureInRegion[x]==0){
        featurelib_length=length(featurelib)+1
        featureInRegion[x]=featurelib_length
        featurelib[[featurelib_length]]<<-feature
      }
    }
  }
  return(featureInRegion)
}

MatrixScan=function(regionMatrix,wd=3,distanceThreshold=0.2){
  #regionMatrix: each row is a position on genome, each col is a marker
  #wd: the length of feature window
  #featurelib is a list of feature that showed up before, must exist in parent enviornment
  wdm=wd-1
  featureInRegion=rep(0,(nrow(regionMatrix)-wdm))
  #the first region
  feature=regionMatrix[1:(1+wdm),]
  featureInRegion[1]=FindFeature(feature,featurelib,distanceThreshold)
  #the following region
  for(x in 2:(nrow(regionMatrix)-wdm)){
    lastfeature=regionMatrix[(x-1):(x+wdm-1),]
    feature=regionMatrix[x:(x+wdm),]
    if(sum(abs(lastfeature-feature))<distanceThreshold){featureInRegion[x]=featureInRegion[x-1]}
    else{featureInRegion[x]=FindFeature(feature,featurelib,distanceThreshold)}
  }
  return(featureInRegion)
}

FindFeature=function(feature,featurelib,distanceThreshold){
  if(length(featurelib)>0){
    distance=sapply(1:length(featurelib),function(x){sum(abs(featurelib[[x]]-feature))})
    tagFeature=which.min(distance)
    if(distance[tagFeature]<distanceThreshold){return(tagFeature)}
    else{return(0)}
  }
  else{return(0)}
}

FeatureDistributionListGenrator_abs=function(positive_sample_list,negative_sample_list,cluster_idx){
  featureDistriList=c()
  posifeatureDistriList=vector(mode="list", length=max(cluster_idx))#initialize list of vector
  posifeaturePresentNum=rep(0,max(cluster_idx))
  negfeatureDistriList=vector(mode="list", length=max(cluster_idx))
  negfeaturePresentNum=rep(0,max(cluster_idx))
  #gather postitional infor
  for(i in 1:length(positive_sample_list)){
    currentFeatureVec=cluster_idx[positive_sample_list[[i]]]
    featurePresent=unique(currentFeatureVec)
    posifeaturePresentNum[featurePresent]=posifeaturePresentNum[featurePresent]+1
    for(k in featurePresent){
      featurePosi=pmin(which(currentFeatureVec==k),length(currentFeatureVec)+1-which(currentFeatureVec==k))
      posifeatureDistriList[[k]]=c(posifeatureDistriList[[k]],featurePosi)
    }
  }
  for(i in 1:length(negative_sample_list)){
    currentFeatureVec=cluster_idx[negative_sample_list[[i]]]
    featurePresent=unique(currentFeatureVec)
    negfeaturePresentNum[featurePresent]=negfeaturePresentNum[featurePresent]+1
    for(k in featurePresent){
      featurePosi=pmin(which(currentFeatureVec==k),length(currentFeatureVec)+1-which(currentFeatureVec==k))
      negfeatureDistriList[[k]]=c(negfeatureDistriList[[k]],featurePosi)
    }
  }
  #stat & smooth
  #two vector not same length
  for(k in 1:max(cluster_idx)){
    vec1=Feature_posistat_vec(posifeatureDistriList[[k]])/posifeaturePresentNum[k]
    vec2=Feature_posistat_vec(negfeatureDistriList[[k]])/negfeaturePresentNum[k]
    if(length(vec1)>length(vec2)){vec2=c(vec2,rep(0,length(vec1)-length(vec2)))}
    if(length(vec1)<length(vec2)){vec1=c(vec1,rep(0,length(vec2)-length(vec1)))}
    featureDistriList[[k]]=vec1-vec2
  }
  return(featureDistriList)
}


#####################################################
#DMD
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



###################################################3
#DTW
shrinked_matrixlist_500=lapply(1:nrow(sampleset_500),function(x){return(matrix_shrink(matrixlist_500[[x]]))})

#distance mastrix generation
dis_matrix_500=dtw_distMatCreat(matrixlist_500)
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/dis_mastrix_50")
pdf(file="test.pdf",width = 15,height = 15)
breaks=seq(0,3,0.1)
heatmap.2(log10(dis_matrix_500+1),col=colorRampPalette(c("yellow","red","black")),
          trace="none",density.info="histogram",Colv=F,Rowv=F,notecol="black",dendrogram = "none",
          lhei=c(1,5),lwid=c(1,5),margins = c(5.5,5.5))
dev.off()

pdf(file="test.pdf",width = 5,height = 5)
rtsne_out2=Rtsne(dis_matrix_500,dims=2,is_distance=T,perplexity=5)
plot(rtsne_out2$Y,col=c(rep(1,25),rep(2,25)))
dev.off()
#tsne mapping

#time testing
ptm <- proc.time()
dtw(matrix_shrink(matrixlist_500[[1]]),matrix_shrink(matrixlist_500[[2]]))$distance
proc.time() - ptm

ptm <- proc.time()
dtw(matrixlist_500[[1]],matrixlist_500[[2]])$distance
proc.time() - ptm



feature build_scan_analysis
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


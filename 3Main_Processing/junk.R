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
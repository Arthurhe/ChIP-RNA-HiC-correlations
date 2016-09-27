#these functions require the following libraries
#library(proxy) #dist
#library(fastcluster)
#for(i in list.files(path="F:/DATA/R/my_functions",full.names = T)){source(i)}

#################################
#function in db_build
get_MarkMatrix=function(position_info,dbhandle,target_cell_type,target_markers){
  tblname_list=paste(target_cell_type,target_markers,sep="_")
  chr=position_info[1]
  posstart=position_info[2]-500
  posstopstart=position_info[3]-500
  coverage_list=lapply(tblname_list,function(x){return(
    dbGetQuery(dbhandle,paste("SELECT coverage FROM ",x," WHERE chr='",chr,"' AND start>",posstart," AND start<",posstopstart,sep="")))})
  marker_matrix=do.call(cbind, coverage_list)
  colnames(marker_matrix)=target_markers
  return(marker_matrix)
}

#still working on the batch stuff
get_MarkMatrix_batch=function(position_info,dbhandle,target_cell_type,target_markers){
  tblname_list=paste(target_cell_type,target_markers,sep="_")
  chr=position_info[,1]
  posstart=position_info[,2]+5000-500
  posstopstart=position_info[,3]+20000-500
  #generate the commande
  for(present_chr in unique(chr)){
    command=paste(" WHERE chr=",present_chr)
    for(numx in 1:nrow(position_info)){
      command=paste(command," chr=",sep="")
    }
    coverage_list=lapply(tblname_list,function(x){return(
      dbGetQuery(dbhandle,paste("SELECT coverage FROM ",x," WHERE chr='",chr,"' AND start>",posstart," AND start<",posstopstart,sep="")))})
  }
  marker_matrix=do.call(cbind, coverage_list)
  rownames(marker_matrix)=target_markers
  return(marker_matrix)
}

matrix_shrink=function(thematrix,threshold=0.1,getridofmid=T){
  if(getridofmid){
    thematrix=thematrix[seq(1,nrow(thematrix),2),]
  }
  therows=rowSums(thematrix)
  wanted=sapply(1:(length(therows)-1),function(x){return(sum(therows[x:(x+1)]))})
  wanted=c(wanted,1)
  return(thematrix[wanted>threshold,])
}

dtw_1toMany=function(theone,themany,distmethod="Manhattan"){
  return(sapply(themany,function(i){dtw(theone,i,dist.method = distmethod)$distance}))
}

dtw_distMatCreat=function(thelist,distmethod="Manhattan")
{
  m = diag(0, length(thelist))
  sapply(1:(length(thelist)-1), function(i)
  {
    m[,i] <<- c(rep(0,i), dtw_1toMany(thelist[[i]],thelist[(i+1):length(thelist)],distmethod))
  }) 
  m= m + t(m)
  return(m)
}

#################################
#function in feature extraction

list <- structure(NA,class="result")
#redifine list so that list can contain empty value
#its espcially useful in function returning multiple value
#check http://stackoverflow.com/questions/1826519/function-returning-more-than-one-value
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

KmeanFeatureClustering=function(featurelib_matrix,centers,nstart,iter.max){
  cluster_info=kmeans(featurelib_matrix,centers,nstart = nstart,iter.max=iter.max)
  truedis=dist(featurelib_matrix,cluster_info$centers,method="manhattan")
  cluster_info$cluster=apply(truedis,1,which.min)
  return(cluster_info)
}

MatrixScan_featureLibBuild_Advance=function(regionMatrix,wd=1,featurelib_matrix=NA,cluster_info=list("cluster"=0),distanceThreshold=0.2){
  #regionMatrix: each row is a position on genome, each col is a marker
  #featurelib_matrix is a matrix feature that showed up before
  #cluster_info is the kmean object returned by kmean featurelib_matrix
  #wd: the length of feature window
  #wd must be identical with feature length in library
  if(is.data.frame(featurelib_matrix) && ncol(featurelib_matrix)/71!=wd){stop('incorrect featurelib_matrix col num')}
  wdm=wd-1
  #craete rawfeature list(in matrix format) for regionMatrix
  if(wd==1){feature=regionMatrix
  }else{
    feature_temp=lapply(1:wd,function(x){regionMatrix[seq(x,nrow(regionMatrix)-wd+x),]})
    feature=do.call(cbind,feature_temp)
  }
  #clustering the new raw feature
  d=dist(feature,method = "manhattan",upper = T) #euclidean
  raw_feature_cluster=hclust(d,method="complete") #complete should be used here
  cluster_idx=as.vector(cutree(raw_feature_cluster, h=distanceThreshold))
  unique_feature=as.data.frame(GetGroupCenter(feature,cluster_idx))
  #find what are the feature that correspond to raw unique_feature
  unique_feature_assignment=FindFeature_Advance2(unique_feature,featurelib_matrix,cluster_info,distanceThreshold)  
  #deal with new feature
  if(min(unique_feature_assignment)==0){
    #anex feature to library
    newfeature_idx=which(unique_feature_assignment==0)
    newfeature=unique_feature[newfeature_idx,]
    if(is.data.frame(featurelib_matrix)){
      feature_starting_point=nrow(featurelib_matrix)
      featurelib_matrix=rbind(featurelib_matrix,newfeature)
    }
    else{
      feature_starting_point=0
      featurelib_matrix=newfeature
    } #if lib is empty
    #update library clustering (index)
    if(nrow(featurelib_matrix)>500){
      if(max(cluster_info$cluster)==as.integer(nrow(featurelib_matrix)/250)){
        cluster_info=KmeanFeatureClustering(featurelib_matrix,centers=cluster_info$centers,nstart=1,iter.max = 10)
      }else{
        cluster_info=KmeanFeatureClustering(featurelib_matrix,centers=as.integer(nrow(featurelib_matrix)/250),nstart=10,iter.max=50)
      }
    }
    #idx new unique_feature
    unique_feature_assignment[newfeature_idx]=seq((feature_starting_point+1),(feature_starting_point+length(newfeature_idx)))
  }
  #idx the original feature list
  featureInRegion=unique_feature_assignment[cluster_idx]
  return(list(featureInRegion,featurelib_matrix,cluster_info))
}

GetGroupCenter=function(sampleMatrix,cluster_idx){
  center_list=lapply(1:max(cluster_idx),function(x){
    tagRows=which(cluster_idx==x)
    if(length(tagRows)>1){
      return(colMeans(sampleMatrix[tagRows,]))
    }else{
      return(sampleMatrix[tagRows,])
    }
  })
  center_matrix=do.call(rbind,center_list)
  return(center_matrix)
}

MatrixScan_Advance2=function(regionMatrix,featurelib_matrix,cluster_info,distanceThreshold=0.2){
  #MatrixScan_cluster base on multi to multi ability of "dist"
  #regionMatrix: each row is a position on genome, each col is a marker
  #featurelib_matrix is a matrix feature that showed up before
  #cluster_info is the kmean object returned by kmean featurelib_matrix
  wd=ncol(featurelib_matrix)/71 #wd: the length of feature window
  wdm=wd-1
  region_length=nrow(regionMatrix)-wdm
  if(wd%%1!=0){stop('incorrect featurelib_matrix col num')}
  if(!(length(featurelib_matrix)>1 && length(cluster_info)>1)){stop('featurelib_matrix & cluster_info must exist')}
  #craete rawfeature list(matrix) for regionMatrix
  if(wd==1){feature=regionMatrix
  }else{
    feature_temp=lapply(1:wd,function(x){regionMatrix[seq(x,nrow(regionMatrix)-wd+x),]})
    feature=do.call(cbind,feature_temp)
  }
  #clustering the new raw feature
  d=dist(feature,method = "manhattan",upper = T) #euclidean
  raw_feature_cluster=hclust(d,method="complete") #complete should be used here
  cluster_idx=as.vector(cutree(raw_feature_cluster, h=distanceThreshold))
  unique_feature=as.data.frame(GetGroupCenter(feature,cluster_idx))
  #find what are the feature that correspond to raw unique_feature
  unique_feature_assignment=FindFeature_Advance2(unique_feature,featurelib_matrix,cluster_info,distanceThreshold)  
  #assign the feature idx to raw feature
  featureInRegion=unique_feature_assignment[cluster_idx]
  return(featureInRegion)
}

FindFeature_Advance2=function(feature,featurelib_matrix=NA,cluster_info=NA,distanceThreshold){
  #feature is a data.frame! of multiple feature by row
  if(length(cluster_info)>1){
    distance=dist(feature,cluster_info$centers,method = "manhattan")
    cluster_idx=apply(distance,1,which.min)
    tagFeature_vec=sapply(1:nrow(feature),function(x){
      pick_feature=which(cluster_info$cluster==cluster_idx[x])
      distance=dist(feature[x,],featurelib_matrix[pick_feature,],method = "manhattan")
      if(min(distance)<distanceThreshold){tagFeature=pick_feature[which.min(distance)];return(as.integer(tagFeature))
      }else{return(0)}
    })
    return(tagFeature_vec)
  }
  else if(length(featurelib_matrix)>1){
    distance=dist(feature,featurelib_matrix,method = "manhattan")
    tagFeature_vec=sapply(1:nrow(feature),function(x){
      distance=dist(feature[x,],featurelib_matrix,method = "manhattan")
      if(min(distance)<distanceThreshold){tagFeature=which.min(distance);return(as.integer(tagFeature))
      }else{return(0)}
    })
    return(tagFeature_vec) 
  }
  else{
    return(rep(0,nrow(feature)))
  }
}

FeatureMatrixGenerator=function(feature_vec_list,featureNum,cluster_idx=NULL){
  #generate feature_matrix, which is a matrix showing times of apperance of each feature in each sample
  #feature_vec_list is a list of samples, each sample represented by a vector of features (the output of feature_vec_list)
  #cluster_idx is the index list of cluster of all feature, calculated by cluster algo
  sampleNum=length(feature_vec_list)
  featureMatrix=matrix(0,sampleNum,featureNum)
  for(i in 1:sampleNum){
    feature_count_table=as.data.frame(table(feature_vec_list[i]))
    feature_count_table[,1]=as.numeric(levels(feature_count_table[,1])[feature_count_table[,1]])#convert factor to numeric
    for(x in 1:nrow(feature_count_table)){featureMatrix[i,feature_count_table[x,1]]=feature_count_table[x,2]}
  }
  if(!is.null(cluster_idx)){
    featureMatrix_merge=matrix(0,sampleNum,max(cluster_idx))
    for(i in 1:max(cluster_idx)){
      pick=which(cluster_idx==i)
      if(length(pick)>1){featureMatrix_merge[,i]=rowSums(featureMatrix[,pick])
      }else{featureMatrix_merge[,i]=featureMatrix[,pick]}
    }
    featureMatrix=featureMatrix_merge
  }
  return(featureMatrix)
}

IDF_calculator=function(featureMatrix){
  #count idf (inverse document frequency) for each feature in featureMatrix
  #featureMatrix is a matrix showing times of apperance of each feature in each sample
  idf=log10(colSums(featureMatrix)/colSums(featureMatrix>0))
  return(idf)
}

TF_calculator=function(featureMatrix){
  #count tf (term frequency) for each feature in featureMatrix
  #featureMatrix is a matrix showing times of apperance of each feature in each sample
  tf=t(apply(featureMatrix,1,function(x){x/sum(x)}))
  return(tf)
}

TF_IDF_calculator=function(tf,idf){
  #count tf-idf (term frequency-inverse document frequency) for each feature in featureMatrix
  tf_idf=t(apply(tf,1,function(x){x*idf}))
  return(tf_idf)
}

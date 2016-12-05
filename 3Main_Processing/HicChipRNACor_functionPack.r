#these functions require the following libraries
#library(proxy) #dist
#library(fastcluster)
#for(i in list.files(path="F:/DATA/R/my_functions",full.names = T)){source(i)}

#################################
#function in db_build
get_MarkMatrix=function(position_info,dbhandle,target_cell_type,target_markers){
  #create the namelist for grepping (cell_type + marker)
  tblname_list=paste(target_cell_type,target_markers,sep="_")
  chr=position_info[1]
  posstart=position_info[2]-500
  posstopstart=position_info[3]
  coverage_list=lapply(tblname_list,function(x){return(
    dbGetQuery(dbhandle,paste("SELECT coverage FROM ",x," WHERE chr='",chr,"' AND start>",posstart," AND start<",posstopstart,sep="")))})
  #the matrix will composed fo marker in coloum, position in row
  marker_matrix=do.call(cbind, coverage_list)
  colnames(marker_matrix)=target_markers
  return(marker_matrix)
}

matrix_percentail_shrink=function(marker_matrix,divide_into=10){
  if(divide_into<1){stop("must divide the matrix into more than 1 chunk")}
  rownum=nrow(marker_matrix)
  merging_pair=matrix(0,divide_into,2)
  merging_pair[,1]=round(seq(1,rownum,rownum/divide_into))
  merging_pair[,2]=round(seq(rownum/divide_into,rownum,rownum/divide_into))
  newmatrix=matrix(0,divide_into,ncol(marker_matrix))
  for(i in 1:divide_into){
    newmatrix[i,]=colMeans(marker_matrix[merging_pair[i,1]:merging_pair[i,2],])
  }
  colnames(newmatrix)=colnames(marker_matrix)
  return(newmatrix)
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

matrix_nozeros_shrink=function(thematrix,threshold=0.1,getridofmid=T){
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

FeatureMatrixGenerator=function(feature_vec_list,featureNum,cluster_idx=1:length(featureNum),peak_region_matrix=NULL){
  #generate feature_matrix, which is a matrix showing times of apperance of each feature in each sample
  #feature_vec_list is a list of samples, each sample represented by a vector of features (the output of feature_vec_list)
  #cluster_idx is the index list of cluster of all feature, calculated by cluster algo
  #region_matrix is a 2 col matrix, representing which region is the region that count for feature appearance
  sampleNum=length(feature_vec_list)
  featureWCenter=as.integer(rownames(peak_region_matrix))
  featureMatrix=matrix(0,sampleNum,max(cluster_idx))
  for(i in 1:sampleNum){
    convert_to_clusteridx=as.vector(cluster_idx[feature_vec_list[[i]]])
    feature_count_table=Vector_stat_table(convert_to_clusteridx)
    rownames(feature_count_table)=feature_count_table[,1]
    featureWCenter_current=feature_count_table[,1][feature_count_table[,1] %in% featureWCenter]
    featureNoCenter_current=feature_count_table[,1][!feature_count_table[,1] %in% featureWCenter]
    for(x in featureWCenter_current){
        featurePosi=which(feature_vec_list[[i]]==x)
        featureMatrix[i,feature_count_table[as.character(x),1]]=sum(featurePosi>=peak_region_matrix[as.character(x),1] & featurePosi<=peak_region_matrix[as.character(x),2])
    }
    for(x in featureNoCenter_current){
        featureMatrix[i,feature_count_table[as.character(x),1]]=feature_count_table[as.character(x),2]
    }
  }
  return(featureMatrix)
}

IDF_calculator=function(featureMatrix,ignore=NULL){
  #count idf (inverse document frequency) for each feature in featureMatrix
  #featureMatrix is a matrix showing times of apperance of each feature in each sample
  idf=log10(colSums(featureMatrix)/colSums(featureMatrix>0))
  idf[ignore]=0
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

FeatureDistributionListGenrator_relative=function(positive_sample_list,cluster_idx){
  featureDistriMatrix=matrix(0,max(cluster_idx),100)
}

FeatureDistriListGenrator_abs=function(sample_list,cluster_idx,Poisson_threshold=0.05){
  featureDistriList=vector(mode="list", length=max(cluster_idx))#initialize list of vector
  for(i in 1:length(sample_list)){
    currentFeatureVec=cluster_idx[sample_list[[i]]]
    featurePresent=unique(currentFeatureVec)
    for(k in featurePresent){
      featurePosi=pmin(which(currentFeatureVec==k),length(currentFeatureVec)+1-which(currentFeatureVec==k)) #the position feature showed up
      featureDistriList[[k]]=c(featureDistriList[[k]],featurePosi)
    }
  }
  #convert all postion number to stat of number of apperance in each position
  featureDistriList=lapply(1:max(cluster_idx),function(k){Feature_posistat_vec(featureDistriList[[k]])})
  return(featureDistriList)
}

FeaturePositionScoreMatrixGenerator=function(featureDistriList_pos,featureDistriList_neg=NULL,Poisson_threshold=0.01){
  #generate the region that matters for positive feature or negative feature
  if(is.null(featureDistriList_neg)){
    peak_region_list=lapply(1:length(featureDistriList_pos),function(k){
      peak_region=FeaturePeakFinder(featureDistriList_pos[[k]],Poisson_threshold=Poisson_threshold)
      return(peak_region)
    })
    peak_region_list=rlist::list.clean(peak_region_list)
    peak_region_matrix=do.call(rbind,peak_region_list)
    rownames(peak_region_matrix)=names(peak_region_list)
  }else{
    peak_region_list=lapply(1:length(featureDistriList_pos),function(k){
      pos_length=length(featureDistriList_pos[[k]])
      neg_length=length(featureDistriList_neg[[k]])
      if(pos_length>neg_length){
        featureDistriVec=featureDistriList_pos[[k]]-c(featureDistriList_neg[[k]],rep(0,pos_length-neg_length))
      }else if(pos_length<neg_length){
        featureDistriVec=c(featureDistriList_pos[[k]],rep(0,neg_length-pos_length))-featureDistriList_neg[[k]]
      }else{featureDistriVec=featureDistriList_pos[[k]]-featureDistriList_neg[[k]]}
      pos_peak_region=FeaturePeakFinder(featureDistriVec,Poisson_threshold=Poisson_threshold)
      return(pos_peak_region)
    })
    names(peak_region_list)=1:length(peak_region_list)
    peak_region_list=rlist::list.clean(peak_region_list)
    peak_region_matrix=do.call(rbind,peak_region_list)
    rownames(peak_region_matrix)=names(peak_region_list)
  }
  return(peak_region_matrix)
}

FeaturePeakFinder=function(featurePosiVec,Poisson_threshold=0.05){
  #find peak in position stat vector, only get positive peak or negative peak (not both) for now
  if(mean(featurePosiVec)<0){featurePosiVec=-featurePosiVec}
  if(qpois(Poisson_threshold,mean(featurePosiVec),lower.tail = F)<=max(featurePosiVec)){#whether there is a peak
    #get consecutive region around peak
    peakPosi=which(featurePosiVec==max(featurePosiVec)) #peak position
    candidate_positions=featurePosiVec>mean(featurePosiVec) #region candidate position that larger than mean
    derivity=candidate_positions-c(0,candidate_positions[-length(candidate_positions)])#the change trend of each position relative to the previous one
    starts=which(derivity>0)-1 #region starts, elasticity 1(1kb)
    ends=which(derivity<0) #region ends
    if(length(ends)!=length(starts)){ends=c(ends,length(featurePosiVec))}
    regions_length_order=order(ends-starts,decreasing = T)
    for(i in regions_length_order){if(any(peakPosi>=starts[i] & peakPosi<=ends[i])){peak_region=c(starts[i],ends[i]); break}}
    return(peak_region)
  }else{
    return(NULL)
  }
}

#not useing now
FeaturePositionScoreCalculator=function(learning_feature_vec,cluster_idx,featureDistriList){
  sampleNum=length(learning_feature_vec)
  featureMatrix=matrix(0,sampleNum,max(cluster_idx))
  for(i in 1:length(learning_feature_vec)){
    currentFeatureVec=cluster_idx[learning_feature_vec[[i]]]
    featurePresent=unique(currentFeatureVec)
    for(k in featurePresent){
      featurePosi=pmin(which(currentFeatureVec==k),length(currentFeatureVec)+1-which(currentFeatureVec==k))
      featureMatrix[i,k]=sum(featureDistriList[[k]][featurePosi],na.rm = T)
    }
  }
  return(featureMatrix)
}

Feature_numstat_vec=function(featurevec,featurenum){
  #return a vector for the stat of features in a feature vector
  feature_tbl=Vector_stat_table(featurevec)
  feature_stat=rep(0,featurenum+1)
  for(x in 2:nrow(feature_tbl)){feature_stat[feature_tbl[x,1]]=feature_tbl[x,2]}
  feature_stat[featurenum+1]=feature_tbl[1,1]
  return(feature_stat)
}

Feature_posistat_vec=function(featurePosiVec){
  if(length(featurePosiVec)==0){return(0)}
  max_length=max(featurePosiVec)
  feature_posi_tbl=Vector_stat_table(featurePosiVec)
  feature_posi_stat=rep(0,max_length)
  for(x in 1:nrow(feature_posi_tbl)){feature_posi_stat[feature_posi_tbl[x,1]]=feature_posi_tbl[x,2]}
  return(feature_posi_stat)
}

######################
#General function 
General_matrix2table=function(target){
  #a general function to convert matrix into table,
  #with rownames of the matrix in the 1st col of table
  #colnames of the matrix in the 2nd col of table
  #values in the 3rd col
  result=lapply(1:nrow(target),function(r){return(cbind(rep(rownames(target)[r],ncol(target)),colnames(target),t(target[r,])))})
  result=do.call(rbind,result)
  result=as.data.frame(result,stringsAsFactors = F)
  result[,3]=as.numeric(result[,3])
  colnames(result)=c("rname","cname","val")
  rownames(result)=NULL
  return(result)
}

substrRight =function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

Vector_stat_table=function(vec){
  vec_tbl=as.data.frame(table(vec))
  vec_tbl[,1] <- as.numeric(as.character(vec_tbl[,1]))
  return(vec_tbl)
}
#Main Analysis 
library(RSQLite)
library(data.table)
#library(gplots)
#library(dtw)
#library(Rtsne) ##Rtsne 
library(proxy)
library(fastcluster)
source("/home/ahe/Analysis/201608_HicChipRnaCor/codes/3Main_Processing")
#source("F:/DATA/R/Kees/1608_HicChipRNACor/3Main_Processing/HicChipRNACor_functionPack.r")
setwd("/home/ahe/Analysis/201608_HicChipRnaCor/data/")

target_cell_type="GM12"
interactionlist_pos=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/GM12878/GM12878_25K_center_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/GM12878/GM12878_25K_negative_loops.txt",stringsAsFactors = F,sep="\t")

target_markers=c("RNA","CAGE","DNase","FAIRE","RRBS","ATF3","BCL3","BCLAF1","BHLHE40","CEBPB","CHD1","CHD2","CREB1","CTCF","CUX1",
                 "E2F4","EGR1","ELF1","ELK1","EP300","ETS1","EZH2","FOS","GABPA","H2AZ","H3K27ac","H3K27me3","H3K36me3","H3K4me1",
                 "H3K4me2","H3K4me3","H3K79me2","H3K9ac","H3K9me3","H4K20me1","JUND","MAFK","MAX","MAZ","MEF2A","MYC","NFE2","NFYA",
                 "NFYB","NR2C2","NRF1","PML","POLR2A","POLR3G","RAD21","RCOR1","REST","RFX5","SIX5","SMC3","SPI1","SP1","SRF","STAT1",
                 "STAT5A","TAF1","TBL1XR1","TBP","USF1","USF2","YY1","ZBTB33","ZNF143","ZNF274","ZNF384","HCFC1")

dbhandle=dbConnect(SQLite(), dbname = '/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/database/tilingdata.sqlite')

#get 1000 samples for feature_lib build
target_pos_list=interactionlist_pos[interactionlist_pos[,4]<200000,]
target_neg_list=interactionlist_neg[interactionlist_neg[,4]<200000,]
pos_samp_num=500
neg_samp_num=500
learning_list_idx=rbind(target_pos_list[sample(nrow(target_pos_list),pos_samp_num),],target_neg_list[sample(nrow(target_neg_list),neg_samp_num),])
learning_list=lapply(1:nrow(learning_list_idx),function(x){return(get_MarkMatrix(learning_list_idx[x,],dbhandle,target_cell_type,target_markers))})
#output matrix samples as csv 
#for(i in 1:length(matrixlist_500)){
#  write.csv(matrixlist_500[[i]],file=paste("sample_matrix_",i,".csv",sep=""))
#}
#historical sampleset: sampleset 50_short_1; sampleset 200_short_1 & 2

#initiate all necessary variable
learning_feature_vec=c()
row_vec=c()
time_vec=c()
featurelib_matrix=NA
cluster_info=list("cluster"=0)
#create feature list
for(i in 1:length(learning_list_idx)){
  ptm <- proc.time()
  list[learning_feature_vec[[i]],featurelib_matrix,cluster_info]=
    MatrixScan_featureLibBuild_Advance(learning_list[[i]],wd=1,featurelib_matrix,cluster_info,distanceThreshold=0.25)
  row_vec[i]=nrow(learning_list[[i]])
  time_vec[i]=(proc.time()-ptm)[3]
  print(paste(i,"|",row_vec[i],"rows |",time_vec[i]))
}

#get 1000 samples for scanning
target_pos_list=interactionlist_pos[interactionlist_pos[,4]<200000,]
target_neg_list=interactionlist_neg[interactionlist_neg[,4]<200000,]
pos_samp_num=500
neg_samp_num=500
scanning_list_idx=rbind(target_pos_list[sample(nrow(target_pos_list),pos_samp_num),],target_neg_list[sample(nrow(target_neg_list),neg_samp_num),])
scanning_list=lapply(1:nrow(scanning_list_idx),function(x){return(get_MarkMatrix(scanning_list_idx[x,],dbhandle,target_cell_type,target_markers))})

#initiate all necessary variable
scanning_feature_vec=c()
scan_row_vec=c()
scan_time_vec=c()
#featurelib_matrix=do.call(rbind,featurelib)
#cluster_info=KmeanFeatureClustering(featurelib_matrix,as.integer(nrow(featurelib_matrix)/250))
for(i in 1:length(scanning_list)){
  ptm <- proc.time()
  scanning_feature_vec[[i]]=MatrixScan_Advance2(scanning_list[[i]],featurelib_matrix,cluster_info,distanceThreshold=0.2)
  scan_row_vec[i]=nrow(scanning_list[[i]])
  scan_time_vec[i]=(proc.time()-ptm)[3]
  print(paste(i,"|",scan_row_vec[i],"rows |",scan_time_vec[i]))
}

#pick target
target_cell_type="K562" 
interactionlist_pos=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/K562/K562_25K_center_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/K562/K562_25K_negative_loops.txt",stringsAsFactors = F,sep="\t")

#get 1000 samples for scanning
target_pos_list=interactionlist_pos[interactionlist_pos[,4]<200000,]
target_neg_list=interactionlist_neg[interactionlist_neg[,4]<200000,]
pos_samp_num=500
neg_samp_num=500
K562_scanning_list_idx=rbind(target_pos_list[sample(nrow(target_pos_list),pos_samp_num),],target_neg_list[sample(nrow(target_neg_list),neg_samp_num),])
K562_scanning_list=lapply(1:nrow(K562_scanning_list_idx),function(x){return(get_MarkMatrix(K562_scanning_list_idx[x,],dbhandle,target_cell_type,target_markers))})

#initiate all necessary variable
K562_scanning_feature_vec=c()
K562_scan_row_vec=c()
K562_scan_time_vec=c()
#featurelib_matrix=do.call(rbind,featurelib)
#cluster_info=KmeanFeatureClustering(featurelib_matrix,as.integer(nrow(featurelib_matrix)/250))
for(i in 1:length(K562_scanning_list_idx)){
  ptm <- proc.time()
  K562_scanning_feature_vec[[i]]=MatrixScan_Advance2(K562_scanning_list[[i]],featurelib_matrix,cluster_info,distanceThreshold=0.2)
  K562_scan_row_vec[i]=nrow(scanning_list[[i]])
  K562_scan_time_vec[i]=(proc.time()-ptm)[3]
  print(paste(i,"|",scan_row_vec[i],"rows |",scan_time_vec[i]))
}

save(featurelib,scanning_feature_vec,K562_scanning_feature_vec,file="featureScanOn200Short2.Rdata")



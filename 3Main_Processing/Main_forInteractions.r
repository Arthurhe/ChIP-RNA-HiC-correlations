#Main Analysis 
library(RSQLite)
library(data.table)
#library(gplots)
#library(dtw)
#library(Rtsne) ##Rtsne 
library(proxy) #dist
library(fastcluster)
source("/home/ahe/Analysis/201608_HicChipRnaCor/codes/3Main_Processing/HicChipRNACor_functionPack.r")
#source("F:/DATA/R/Kees/1608_HicChipRNACor/3Main_Processing/HicChipRNACor_functionPack.r")
setwd("/home/ahe/Analysis/201608_HicChipRnaCor/data/")

target_cell_type="GM12"
interactionlist_pos=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/GM12878/GM12878_25K_center_le200K_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/GM12878/GM12878_25K_nonOverlap_le200K_negative_loops.txt",stringsAsFactors = F,sep="\t") #GM12878_25K_negative_loops.txt

target_markers=c("RNA","CAGE","DNase","FAIRE","RRBS","ATF3","BCL3","BCLAF1","BHLHE40","CEBPB","CHD1","CHD2","CREB1","CTCF","CUX1",
                 "E2F4","EGR1","ELF1","ELK1","EP300","ETS1","EZH2","FOS","GABPA","H2AZ","H3K27ac","H3K27me3","H3K36me3","H3K4me1",
                 "H3K4me2","H3K4me3","H3K79me2","H3K9ac","H3K9me3","H4K20me1","JUND","MAFK","MAX","MAZ","MEF2A","MYC","NFE2","NFYA",
                 "NFYB","NR2C2","NRF1","PML","POLR2A","POLR3G","RAD21","RCOR1","REST","RFX5","SIX5","SMC3","SPI1","SP1","SRF","STAT1",
                 "STAT5A","TAF1","TBL1XR1","TBP","USF1","USF2","YY1","ZBTB33","ZNF143","ZNF274","ZNF384","HCFC1")

dbhandle=dbConnect(SQLite(), dbname = '/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/database/tilingdata.sqlite')

#get 1000 samples for feature_lib build
target_pos_list=interactionlist_pos[interactionlist_pos[,4]<200000,]
target_neg_list=interactionlist_neg[interactionlist_neg[,4]<200000,]
pos_samp_num=5000
neg_samp_num=5000
learning_list_idx=rbind(target_pos_list[sample(nrow(target_pos_list),pos_samp_num),],target_neg_list[sample(nrow(target_neg_list),neg_samp_num),])
learning_list=lapply(1:nrow(learning_list_idx),function(x){return(get_MarkMatrix(learning_list_idx[x,],dbhandle,target_cell_type,target_markers))})

#scan all shit
allGMshort_list_idx=rbind(target_pos_list,target_neg_list)
allGMshort_list=lapply(1:nrow(allGMshort_list_idx),function(x){return(get_MarkMatrix(allGMshort_list_idx[x,],dbhandle,target_cell_type,target_markers))})
allGMshort_list=learning_list
interactionlist_pos=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/GM12878/GM12878_25K_center_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/GM12878/GM12878_25K_negative_loops.txt",stringsAsFactors = F,sep="\t") #GM12878_25K_negative_loops.txt
target_pos_list=interactionlist_pos[interactionlist_pos[,4]<1000000,]
target_neg_list=interactionlist_neg[interactionlist_neg[,4]<1000000,]
pos_samp_num=5000
neg_samp_num=5000
GM1M_list_idx=rbind(target_pos_list[sample(nrow(target_pos_list),pos_samp_num),],target_neg_list[sample(nrow(target_neg_list),neg_samp_num),])
GM1M_list=lapply(1:nrow(GM1M_list_idx),function(x){return(get_MarkMatrix(GM1M_list_idx[x,],dbhandle,target_cell_type,target_markers))})

save(allGMshort_list,GM1M_list,file="matrix_allGM_short.Rdata")


#output matrix samples as csv 
#for(i in 1:length(matrixlist_500)){
#  write.csv(matrixlist_500[[i]],file=paste("sample_matrix_",i,".csv",sep=""))
#}
#historical sampleset: sampleset 50_short_1; sampleset 200_short_1 & 2

#check whether it work on other cell type
target_cell_type="K562" 
interactionlist_pos=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/K562/K562_25K_center_le200K_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/K562/K562_25K_nonOverlap_le200K_negative_loops.txt",stringsAsFactors = F,sep="\t")

#get 5000 samples for scanning
target_pos_list=interactionlist_pos[interactionlist_pos[,4]<200000,]
target_neg_list=interactionlist_neg[interactionlist_neg[,4]<200000,]
pos_samp_num=2500
neg_samp_num=2500
K562_scanning_list_idx=rbind(target_pos_list[sample(nrow(target_pos_list),pos_samp_num),],target_neg_list[sample(nrow(target_neg_list),neg_samp_num),])
K562_scanning_list=lapply(1:nrow(K562_scanning_list_idx),function(x){return(get_MarkMatrix(K562_scanning_list_idx[x,],dbhandle,target_cell_type,target_markers))})

save(learning_list,K562_scanning_list,file="matrix_15000_short.Rdata")

#load("matrix_3000_short.Rdata")
#initiate all necessary variable
learning_feature_vec=c()
scanning_feature_vec=c()
extrascanning_feature_vec=c()
K562_scanning_feature_vec=c()
featurelib_matrix=NA
cluster_info=list("cluster"=0)
#create feature list
for(i in 1:2000){
  list[learning_feature_vec[[i]],featurelib_matrix,cluster_info]=
    MatrixScan_featureLibBuild_Advance(learning_list[[i]],wd=1,featurelib_matrix,cluster_info,distanceThreshold=0.25)
  print(i)
}
for(i in 2001:10000){
  learning_feature_vec[[i]]=MatrixScan_Advance2(learning_list[[i]],featurelib_matrix,cluster_info,distanceThreshold=0.25)
  print(i+2000)
}

for(i in 1:length(K562_scanning_list)){
  K562_scanning_feature_vec[[i]]=MatrixScan_Advance2(K562_scanning_list[[i]],featurelib_matrix,cluster_info,distanceThreshold=0.25)
  print(i+10000)
}
save(featurelib_matrix,cluster_info,learning_feature_vec,K562_scanning_feature_vec,file="featureOn15000Short.Rdata")


#scan all shit
#load("featureOn15000Short.Rdata")
scanning_feature_vec=c()
scanning_long_feature_vec=c()
#create feature list
for(i in 1:length(allGMshort_list)){
  scanning_feature_vec[[i]]=MatrixScan_Advance2(allGMshort_list[[i]],featurelib_matrix,cluster_info,distanceThreshold=0.25)
  print(i)
}

for(i in 1:length(GM1M_list)){
  scanning_long_feature_vec[[i]]=MatrixScan_Advance2(GM1M_list[[i]],featurelib_matrix,cluster_info,distanceThreshold=0.25)
  print(-i)
}
save(featurelib_matrix,cluster_info,scanning_feature_vec,scanning_long_feature_vec,file="featureOnAllGM_short.Rdata")

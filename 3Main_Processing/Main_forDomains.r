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
domain_pos=fread("/home/ahe/Analysis/201608_HicChipRnaCor/data/mydomain/GM_all_elemental_domains.bed",stringsAsFactors = F,sep="\t",select = 1:3,data.table = F)
domain_neg=fread("/home/ahe/Analysis/201608_HicChipRnaCor/data/mydomain/GM_neg_elemental_domains.bed",stringsAsFactors = F,sep="\t",select = 1:3,data.table = F) #GM12878_25K_negative_loops.txt
colnames(domain_neg)=colnames(domain_pos)

target_markers=c("RNA","CAGE","DNase","FAIRE","RRBS","ATF3","BCL3","BCLAF1","BHLHE40","CEBPB","CHD1","CHD2","CREB1","CTCF","CUX1",
                 "E2F4","EGR1","ELF1","ELK1","EP300","ETS1","EZH2","FOS","GABPA","H2AZ","H3K27ac","H3K27me3","H3K36me3","H3K4me1",
                 "H3K4me2","H3K4me3","H3K79me2","H3K9ac","H3K9me3","H4K20me1","JUND","MAFK","MAX","MAZ","MEF2A","MYC","NFE2","NFYA",
                 "NFYB","NR2C2","NRF1","PML","POLR2A","POLR3G","RAD21","RCOR1","REST","RFX5","SIX5","SMC3","SPI1","SP1","SRF","STAT1",
                 "STAT5A","TAF1","TBL1XR1","TBP","USF1","USF2","YY1","ZBTB33","ZNF143","ZNF274","ZNF384","HCFC1")

dbhandle=dbConnect(SQLite(), dbname = '/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/database/tilingdata.sqlite')

#get marker distribution on domains
db_unitsize=1000
db_overlap=500
extension_window_size=1000
extension_window_num=10
devide_into=20
learning_list_idx=rbind(domain_pos,domain_neg)
learning_list=lapply(1:nrow(learning_list_idx),function(x){
  #update idx to extended idx
  grepping_idx=learning_list_idx[x,]
  grepping_idx[2]=grepping_idx[2]-extension_window_size*extension_window_num
  grepping_idx[3]=grepping_idx[3]+extension_window_size*extension_window_num
  #grep raw matrix from db
  raw_matrix=get_MarkMatrix(grepping_idx,dbhandle,target_cell_type,target_markers)
  #extension_window_merging
  rownum=nrow(raw_matrix)
  shrinkrownum=devide_into+extension_window_num*2
  shrinked_matrix=matrix(0,shrinkrownum,ncol(raw_matrix))
  #pre_domain matrix
  for(i in 1:extension_window_num){
    shrinked_matrix[i,]=colSums(raw_matrix[(i*2-1):(i*2+1),])
  }
  extension_window_num_in_raw=extension_window_num*2+1
  shrinked_matrix[(extension_window_num+1):(extension_window_num+devide_into),]=matrix_percentail_shrink(raw_matrix[extension_window_num_in_raw:(rownum-extension_window_num*2),],devide_into)
  for(i in (extension_window_num+devide_into+1):shrinkrownum){
    shrinked_matrix[i,]=colSums(raw_matrix[(rownum-extension_window_num*2+i*2-2):(rownum-extension_window_num*2+i*2),])
  }    
  return(shrinked_matrix)
})
save(learning_list,file="markers_on_domain.Rdata")

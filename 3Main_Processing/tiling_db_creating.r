#create sqlite tables/databases£ºthrough R 
library(RSQLite)
library(data.table)

setwd("/home/ahe/Analysis/201608_HicChipRnaCor/data/")
#create & index
Db_InitNIndex=function(dbhandle,tblname,inputdata,rewrite=T){
  if(dbExistsTable(dbhandle,tblname)){
    dbRemoveTable(dbhandle,tblname);print(paste(tblname,"exist!"))
    }
  dbSendQuery(conn = dbhandle,paste("CREATE TABLE",tblname,"(chr TEXT,start INT,stop INT,count INT,coverage REAL)"))
  if(is.character(inputdata)){
    dbWriteTable(conn=dbhandle, name=tblname, value=inputdata, sep='\t',header=F, append=T)
  }else{
    dbWriteTable(conn=dbhandle, name=tblname, value=inputdata, append=T)
  }
  dbSendQuery(dbhandle,paste("CREATE INDEX ",tblname,"_index on ",tblname," (chr, start)",sep=""))
}

#create target lists
target_cell_type=c("GM12","K562")
target_markers=c("RNA","CAGE","DNase","FAIRE","RRBS","ATF3","BCL3","BCLAF1","BHLHE40","CEBPB","CHD1","CHD2","CREB1","CTCF","CUX1",
                 "E2F4","EGR1","ELF1","ELK1","EP300","ETS1","EZH2","FOS","GABPA","H2AZ","H3K27ac","H3K27me3","H3K36me3","H3K4me1",
                 "H3K4me2","H3K4me3","H3K79me2","H3K9ac","H3K9me3","H4K20me1","JUND","MAFK","MAX","MAZ","MEF2A","MYC","NFE2","NFYA",
                 "NFYB","NR2C2","NRF1","PML","POLR2A","POLR3G","RAD21","RCOR1","REST","RFX5","SIX5","SMC3","SPI1","SP1","SRF","STAT1",
                 "STAT5A","TAF1","TBL1XR1","TBP","USF1","USF2","YY1","ZBTB33","ZNF143","ZNF274","ZNF384","HCFC1")

#create library
dbhandle=dbConnect(SQLite(), dbname = '/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/database/tilingdata.sqlite')
for(i in target_cell_type){
  for(j in target_markers){
    tblname=paste(i,j,sep="_")
    print(tblname)
    #find each target's position
    match_pattern=paste(i,".*",j,"_hg19_1K_tiling.bed",sep="")
    tiling_ChIP = list.files(path="/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/tiling_files", full.names = T, pattern = match_pattern)
    tiling_RNA = list.files(path="/home/ahe/Analysis/201608_HicChipRnaCor/data/RNA/tiling_files", full.names = T, pattern = match_pattern)
    tiling_RRBS = list.files(path="/home/ahe/Analysis/201608_HicChipRnaCor/data/RRBS/tiling_files", full.names = T, pattern = match_pattern)
    tiling_list=c(tiling_ChIP,tiling_RNA,tiling_RRBS)
    if(length(tiling_list)>1){ # if there are more than one candidate then merge
      temp_list=c()
      for(k in 1:length(tiling_list)){
        temp_list[[k]]=fread(tiling_list[k])
      }
      count=temp_list[[1]][,V4]
      coverage=temp_list[[1]][,V5]
      for(k in 2:length(tiling_list)){
        count=count+temp_list[[k]][,V4]
        coverage=coverage+temp_list[[k]][,V5]
      }
      merged_tbl=cbind(temp_list[[1]][,.(V1,V2,V3)],count/length(tiling_list),coverage/length(tiling_list))
      Db_InitNIndex(dbhandle,tblname,merged_tbl)
    }else{ #else direct input
      Db_InitNIndex(dbhandle,tblname,tiling_list[1])
    }
  }  
}


#match_pattern=paste("(",paste0(target_cell_type,collapse="|"),").*(",paste0(target_markers,collapse="|"),")_hg19_1K_tiling.bed",sep="")

#delete index
for(i in "the wannted list"){
  dbSendQuery(db,paste(" DROP INDEX ",tbl_list[[i]],"_index",sep=""))
}

#test timing
library("RSQLite")
db=dbConnect(SQLite(), dbname = '/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/tilingdata.sqlite')
tbl_list=dbListTables(db)
ptm <- proc.time()
temp=c()
for(i in 1:20){
  temp[[i]]=dbGetQuery(db,paste("SELECT * FROM",tbl_list[[i]],"WHERE chr='chr1' AND start>1000000 AND start<4000000"))
}
proc.time() - ptm

#test result
#searching of 3 conditions (chr='chr1' AND start>1000000 AND start<4000000) on 20 list, speed increase 100 fold, from 10s to 0.1s. file size increase from 2.7GB to 4GB. (by indexing on two column)
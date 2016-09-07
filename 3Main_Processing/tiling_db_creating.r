#create sqlite tables/databases£ºthrough R 
#create target lists
library("RSQLite")
target_cell_type=c("GM12","K562")
target_markers=c("RNA","CAGE","DNase","FAIRE","RRBS","H2AZ","H3K27me3","H3K36me3","H3K4me1","H3K4me2","H3K4me3","H3K79me2","H3K9ac","H3K9me3","H4K20me1",
                 "EZH2","POLR2A","REST","CTCF","RAD21","CUX1","ZNF384","SPI1","SP1","RCOR1","MAX","HCFC1","CEBPB","JUND","TBP","SRF")
match_pattern=paste("(",paste0(target_cell_type,collapse="|"),").*(",paste0(target_markers,collapse="|"),")_hg19_1K_tiling.bed",sep="")
tiling_ChIP = list.files(path="/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/tiling_files", full.names = T, pattern = match_pattern)
tiling_RNA = list.files(path="/home/ahe/Analysis/201608_HicChipRnaCor/data/RNA/tiling_files", full.names = T, pattern = match_pattern)
tiling_RRBS = list.files(path="/home/ahe/Analysis/201608_HicChipRnaCor/data/RRBS/tiling_files", full.names = T, pattern = match_pattern)
tiling_list=c(tiling_ChIP,tiling_RNA,tiling_RRBS)

#create library
db=dbConnect(SQLite(), dbname = '/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/database/tilingdata.sqlite')
for(i in tiling_list){
  filename=basename(file_path_sans_ext(i))
  if(dbExistsTable(db,filename)){dbRemoveTable(db,filename)}
  dbSendQuery(conn = db,paste("CREATE TABLE",filename,"
                              (chr TEXT,
                              start INT,
                              stop INT,
                              count INT,
                              coverage REAL)"))
  dbWriteTable(conn=db, name=filename, value=i, sep='\t',header=F, append=T)
}

#create index
for(i in 1:20){
  temp[[i]]=dbSendQuery(db,paste("
                                 CREATE INDEX ",tbl_list[[i]],"_index on ",tbl_list[[i]]," (chr, start)",sep=""))
}

#delete index
for(i in 1:20){
  temp[[i]]=dbSendQuery(db,paste("
                                 DROP INDEX ",tbl_list[[i]],"_index",sep=""))
}

#test timing
library("RSQLite")
db=dbConnect(SQLite(), dbname = '/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/tilingdata.sqlite')
tbl_list=dbListTables(db)
ptm <- proc.time()
temp=c()
for(i in 1:20){
  temp[[i]]=dbGetQuery(db,paste("
                                SELECT coverage
                                FROM",tbl_list[[i]],
                                "WHERE chr='chr1' AND start>1000000 AND start<4000000"))
}
proc.time() - ptm

#test result
#searching of 3 conditions (chr='chr1' AND start>1000000 AND start<4000000) on 20 list, speed increase 100 fold, from 10s to 0.1s. file size increase from 2.7GB to 4GB. (by indexing on two column)
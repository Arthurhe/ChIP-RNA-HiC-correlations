#Main Analysis 
library(RSQLite)
library(data.table)
library(dtw)
library(Rtsne) ##Rtsne 
setwd("/home/ahe/Analysis/201608_HicChipRnaCor/data/")

#functions
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
  posstart=position_info[,2]-500
  posstopstart=position_info[,3]-500
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
  return(sapply(themany,function(i){dtw(matrix_shrink(theone),matrix_shrink(i),dist.method = distmethod)$distance}))
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

#pick target
target_cell_type="K562"
interactionlist_pos=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/K562/K562_25K_center_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/K562/K562_25K_negative_loops.txt",stringsAsFactors = F,sep="\t")

target_cell_type="GM12"
interactionlist_pos=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/GM12878/GM12878_25K_center_interactions_simplified.txt",stringsAsFactors = F,sep="\t")
interactionlist_neg=read.table("/home/ahe/Analysis/201608_HicChipRnaCor/data/HiC/GM12878/GM12878_25K_negative_loops.txt",stringsAsFactors = F,sep="\t")

target_markers=c("RNA","CAGE","DNase","FAIRE","RRBS","ATF3","BCL3","BCLAF1","BHLHE40","CEBPB","CHD1","CHD2","CREB1","CTCF","CUX1",
                 "E2F4","EGR1","ELF1","ELK1","EP300","ETS1","EZH2","FOS","GABPA","H2AZ","H3K27ac","H3K27me3","H3K36me3","H3K4me1",
                 "H3K4me2","H3K4me3","H3K79me2","H3K9ac","H3K9me3","H4K20me1","JUND","MAFK","MAX","MAZ","MEF2A","MYC","NFE2","NFYA",
                 "NFYB","NR2C2","NRF1","PML","POLR2A","POLR3G","RAD21","RCOR1","REST","RFX5","SIX5","SMC3","SPI1","SP1","SRF","STAT1",
                 "STAT5A","TAF1","TBL1XR1","TBP","USF1","USF2","YY1","ZBTB33","ZNF143","ZNF274","ZNF384","HCFC1")

dbhandle=dbConnect(SQLite(), dbname = '/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/database/tilingdata.sqlite')

#get 500 samples to get a taste of dtw distribution in space
#generate name list
pos_samp_num=25
neg_samp_num=25
sampleset_500=rbind(interactionlist_pos[sample(nrow(interactionlist_pos),pos_samp_num),],interactionlist_neg[sample(nrow(interactionlist_neg),neg_samp_num),])
matrixlist_500=lapply(1:nrow(sampleset_500),function(x){return(get_MarkMatrix(sampleset_500[x,],dbhandle,target_cell_type,target_markers))})
shrinked_matrixlist_500=lapply(1:nrow(sampleset_500),function(x){return(matrix_shrink(matrixlist_500[[x]]))})

#distance mastrix generation
dis_matrix_500=dtw_distMatCreat(matrixlist_500)
load("F:/DATA/R/Kees/1608_HicChipRNACor/data/dis_mastrix_50")
breaks=seq(4,5,0.1)
heatmap.2(log10(dis_matrix_500+1),col=colorRampPalette(c("yellow","red","black")),
          trace="none",density.info="histogram",Colv=F,Rowv=F,notecol="black",dendrogram = "none",
          lhei=c(1,5),lwid=c(1,5),margins = c(5.5,5.5))

#tsne mapping
rtsne_out2=Rtsne(dis_matrix_500,dims=2,is_distance=T,perplexity=15)
plot(rtsne_out2$Y,col=c(rep(1,25),rep(2,25)))

#time testing
ptm <- proc.time()
dtw(matrix_shrink(matrixlist_500[[1]]),matrix_shrink(matrixlist_500[[2]]))$distance
proc.time() - ptm

ptm <- proc.time()
dtw(matrixlist_500[[1]],matrixlist_500[[2]])$distance
proc.time() - ptm

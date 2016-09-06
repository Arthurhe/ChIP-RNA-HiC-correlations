library("RSQLite")
library("tools")
db=dbConnect(SQLite(), dbname = '/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/tilingdata.sqlite')
tiling_list=list.files(path="/home/ahe/Analysis/201608_HicChipRnaCor/data/ChIPnlike/tiling_files",full.names = T,pattern = "(K562|GM12).*_hg19_1K_tiling.bed")
record_num=0
for(i in tiling_list){
	record_num=record_num+1
	filename=basename(file_path_sans_ext(i))
	if(dbExistsTable(db,filename)){dbRemoveTable(db,filename)}
	print(paste(record_num,i,sep="    "))
	dbSendQuery(conn = db,paste("CREATE TABLE",filename,"(chr TEXT,start INT,stop INT,count INT,coverage REAL)"))
	dbWriteTable(conn=db, name=filename, value=i, sep='\t',header=F, append=T)
}

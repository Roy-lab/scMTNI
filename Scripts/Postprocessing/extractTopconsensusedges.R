#!/user/bin/env Rscript
##This is Shilu's script with SR updated to take an output path 
#install.packages("stringr", repos='http://cran.us.r-project.org')
library(stringr)
args=commandArgs(trailingOnly = T);
if(length(args)==0)
{
  stop("At least one argument must be supplied(inputfile).n",call. = FALSE)
}else if(length(args)>0)
{
  print("Start!")
}

library(data.table)
#cl=args[1] ## cl="C1"
#n=args[2] ## n=cf0.5 or top5000
#path=args[2] 
cf=args[1] #"top" or "cf"
n=as.numeric(args[2])
outpath=args[3]
    

cellfile=("/mnt/dv/wid/projects7/Roy-singlecell2/ashtonlab/results/H72_de/scmtni/merlin/input/celltype_order.txt")
cellsf=read.table(cellfile,stringsAsFactors = F)
cells=cellsf$V1

data=NULL
for(cell in cells){
  print(cell)
  path=paste0("/mnt/dv/wid/projects7/Roy-singlecell2/ashtonlab/results/H72_de/scmtni/merlin/subsample/analysis/")
  d=fread(paste0(path,cell,"/consensus_edges.txt"))
  d=d[order(d$V3,decreasing = T),]
  if(cf=="top"){
    d1=d[1:n,]
  }else{
    id=which(d$V3>=n)
    d1=d[id,]
  }
  print(paste(cell,nrow(d1)))
  fwrite(d1,file=paste0(outpath,"/",cell,"/consensus_edges_",cf,n,".txt"),quote=F,sep="\t",row.names = F,col.names = F)
  #d1=fread(paste0(path,cell,"/consensus_edges_top",n,".txt"))
  print(paste0(path,"/",cell,"/consensus_edges_",cf,n,".txt"))
  #cellclust=paste0(cluster,"_",cell)
  names(d1)=c("regulator","target",cell)
  if(is.null(data)){
    data=d1
  }else{
    data=merge(data,d1,by=c("regulator","target"),all=T)
  }
}
data[is.na(data)]=0
data$edge=paste0(data$regulator,"->",data$target)
#cells2$V3 <- str_c(cells2$V1,"_", cells2$V2)
#cname=cells2$V3
#selcols=c("edge",cname)
selcols=c("edge",cells)
dall=data[,..selcols]
fwrite(dall,file=paste0(outpath,"/consensus_edges_merged_",cf,n,".txt"),quote=F,sep="\t",row.names = F,col.names = T)



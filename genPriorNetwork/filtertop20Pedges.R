#!/user/bin/env Rscript
## June 15, 2021: filter prior network top get top 20%
args=commandArgs(trailingOnly = T);
if(length(args)==0)
{
  stop("At least one argument must be supplied(inputfile).n",call. = FALSE)
}else if(length(args)>0)
{
  print("Start making plot of comparing predictons to hic data.")
}
indir0=args[1]
indir=args[2]
library(data.table)
#indir0="/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/data/liger_sqrt_ncell50_k8_filterhumanbc_FBS/"
cells=read.table(paste0(indir0,"/celltype_order.txt"),stringsAsFactors = F)
cells=cells$V1
if (substr(indir,nchar(indir),nchar(indir))=="/")
{
  outdir=paste0(substr(indir,1,nchar(indir)-1),"_top0.2/")
}else{
  outdir=paste0(indir,"_top0.2/")
}
system(paste0("mkdir -p ",outdir))
dsum=NULL
library(dplyr)
for(cell in cells){
  d=fread(paste0(indir,cell,"_network.txt"))
  print(paste0(indir,cell,"_network.txt"))
  tfs=unique(d$V1)
  genes=unique(d$V2)
  n=round(0.2*length(tfs)*length(genes))
  d=d[order(d$V3,decreasing = T),]
  if(n<nrow(d)){
    d1=d[1:n,]
  }else{
    d1=d
    print(paste(cell,n,"<",nrow(d)))
    n=nrow(d)
  }
  print(paste(cell,n))
  fwrite(d1,file=paste0(outdir,cell,"_network.txt"),quote=F,sep="\t",row.names = F,col.names = F)
  tfs1=unique(d1$V1)
  genes1=unique(d1$V2)
  dd=data.frame(cell,total=nrow(d),ntf=length(tfs),ngene=length(genes),n,ntf1=length(tfs1),ngene1=length(genes1))
  dreg=d1%>%group_by(V1) %>% summarise(ntarget=n())
  print(summary(dreg$ntarget))
  dsum=rbind(dsum,dd)
}
print(dsum)
write.csv(dsum,file=paste0(outdir,"summary_priornetwork.csv"),row.names=F,quote=F)


library(data.table)
library(dplyr)
source("Scripts/Network_Analysis/genNetall_multi.R")




cvnpath="Results_subsample/analysis/"



cells=read.table("ExampleData/celltype_order.txt",stringsAsFactors = F)

args = commandArgs(trailingOnly=TRUE)
k=args[1]


cells=cells$V1
outdir=paste0(cvnpath,"lda_TFcellbygene/k",k,"/")
print("Outdir is")
print(outdir)

filter="_cf0.8"
## get most divergent topic:
dstat=NULL
for(tp in 1:k){
  for(cell in cells){
    tfs=read.table(paste0(outdir,"TFs_topicid_",cell,".txt"))
    tf=tfs[which(tfs$V2==tp),]
    prefix=paste0("lda_topic",tp,"_",cell)
    genes=read.table(paste0(outdir,"genes_topicid.txt"))
    gene=genes[which(genes$V2==tp),]
    dstat=rbind(dstat,data.frame(tp,cell,tf=nrow(tf),gene=nrow(gene)))
  }
}


dt=dstat %>% group_by(tp) %>% summarise(mean=mean(tf),variance=var(tf)) %>% arrange(desc(variance))
write.csv(dt,file=paste0(outdir,"lda_k",k,"_variance.csv"),quote=F,row.names=F)

## plot networks for each topic:
genes=read.table(paste0(outdir,"genes_topicid.txt"))
names(genes)=c("target","targetid","targetprob")
for(tp in 1:k){ #c(2,3,10)
  ## extract networks for regulators and targets for each topic:
  dtopic=NULL
  gene=genes[which(genes$targetid==tp),]
  for(cell in cells){
    tfs=read.table(paste0(outdir,"TFs_topicid_",cell,".txt"))
    names(tfs)=c("regulator","tfid","tfprob")
    tf=tfs[which(tfs$tfid==tp),]
    print(paste(cell,"TF:",unique(tf$tfid),"gene:",unique(gene$targetid)))
    prefix=paste0("lda_topic",tp,"_",cell)
    d=fread(paste0(cvnpath,"/",cell,"/consensus_edges",filter,".txt"))
    names(d)=c("regulator","target","weight")
    print(nrow(d))
    d1=merge(d,tfs,by="regulator")
    d1=merge(d1,genes,by="target")
    #if(nrow(d)!=nrow(d1)){
      #print("erorr!")
    #}
    #id=which(d1$regulator%in%tf$regulator & d1$target%in%gene$target)
    id=which(d1$tfid==tp & d1$targetid==tp)
    d2=d1[id,]
    d2$cell=cell
    dtopic=rbind(dtopic,d2)
    
  }
  dtopic1=dtopic[,c("regulator","target","weight","cell")]
  fwrite(dtopic1,file=paste0(outdir,"networks_lda_k",k,"_topic",tp,".txt"),quote=F,row.names=F,col.names=T,sep="\t")
  prefix=paste0("lda_topic",tp)
  outname=paste0(outdir,prefix,"_consensus_edges_k",k,"_network.png")
  genNetall_multi(prefix,dtopic1,outname,w=5000,h=3000,nr=2,ntop=20, cells)
}



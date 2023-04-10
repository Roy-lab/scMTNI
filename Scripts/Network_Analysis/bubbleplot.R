library(tidyr)
library(data.table)
library(dplyr)
library(gridExtra)
library(ggplot2)



# Command line arguments
# 1 - number of clusters
# 2 - path to clusters
# 3 - path to outfile
# 4 - prefix of clusters


args = commandArgs(trailingOnly=TRUE)


k=args[1]
path0=args[2]
outpath0=args[3]
prefix0=args[4]

n=5000

d=read.delim(paste0(path0,'/',prefix0,".txt"),stringsAsFactors = F)
cells=colnames(d);
cells=cells[2:(length(cells)-1)];
edge=do.call(rbind,strsplit(d$Edge,"->"))
d$regulator=edge[,1]
d$target=edge[,2]
outpath=paste0(outpath0,"/","/kmeans_k",k,"/")
system(paste("mkdir -p",outpath))
#setwd(outpath)

## make heatmap using average profile of kmeans clustering:
dmean=NULL
for(i in 1:k){ #1:k c(5,4,21) c(1,19,15,12)
  data=d[which(d$cluster==i),]
  prefix=paste0("E",i)
  dmean=rbind(dmean,data.frame(prefix,t(colMeans(data[,2:(length(cells)+1)])),nedge=nrow(data)))
}
row.names(dmean)=dmean$prefix
dmean1=dmean[,-1] %>% arrange_all(desc)
level=row.names(dmean1)
dmean1$prefix=row.names(dmean1)
dmean2=dmean1 %>% pivot_longer(1:length(cells),names_to = "cell", values_to = "score")
dmean2$prefix=factor(dmean2$prefix,levels=rev(level))
dmean2$cell=gsub("cluster","C",dmean2$cell)
dmean2$cell <- factor(dmean2$cell, levels = gsub("cluster","C",cells))

gg=ggplot(dmean2, aes(cell, prefix)) + geom_tile(aes(fill = score),colour = "white") + geom_text(aes(label = round(score,digits = 2)))+
  scale_x_discrete(position = "top")+#coord_equal()+ylim(rev(levels(data$cell1)))+
  theme(text = element_text(size=10 ))+
  theme(axis.text.x = element_text(angle=90,vjust = 1.5,hjust = 0))+
  scale_fill_gradient(low = "white",high = "red",na.value = 'black')+theme(legend.position="bottom")
#gg <- gg + facet_wrap(~algorithm, nrow=1)#scales="free_y")
dmean3=dmean1[,c("prefix","nedge")]
dmean3$total="numEdge"
dmean3$prefix=factor(dmean3$prefix,levels=rev(level))
gg2=ggplot(dmean3, aes(total, prefix)) + geom_tile(aes(fill = nedge),colour = "white") + geom_text(aes(label = round(nedge)))+
  scale_x_discrete(position = "top")+#coord_equal()+ylim(rev(levels(data$cell1)))+
  theme(text = element_text(size=10))+
  scale_fill_gradient(low = "white",high = "purple")+theme(legend.position="bottom")

#setwd(path)
library("grid")
pdf(paste0(outpath,"kmeans_meanprofile_heatmap.pdf"),height=10,width=10)
grid.newpage()
print(gg, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(gg2, vp = viewport(x = 0.90, y = 0.5, width = 0.2, height = 1.0))
dev.off()

dsum=NULL
topRegs=NULL
rtop=5
for(i in 1:k){ #c(2,5,10)
  data=d[which(d$cluster==i),]
  prefix=paste0("E",i)
  ddd=data.frame(prefix,nreg=length(unique(data$regulator)),ntarget=length(unique(data$target)),nedge=nrow(data))
  dsum=rbind(dsum,ddd)
  dreg=data%>%group_by(regulator) %>% summarise(ntarget=n(),target=paste(target,collapse = ";")) %>% arrange(desc(ntarget)) %>% as.data.frame()
  
  if(nrow(dreg)<=rtop){
    dreg1=dreg[,1:2]
  }else{
    cf=dreg[rtop,2]
    dreg1=dreg[which(dreg$ntarget>=cf),1:2]  #dreg[1:rtop,1:2]  #dreg[which(dreg$ntarget>=cf),1:2]
  }
  print(paste(prefix,nrow(dreg1)))
  
  names(dreg1)[2]=prefix
  if(is.null(topRegs)){
    topRegs=dreg1[,1:2]
  }else{
    topRegs=merge(topRegs,dreg1[,1:2],by="regulator",all=T)
  }
}
rmax=apply(topRegs[,-1],1,function(x) max(x,na.rm = T))
rid=which(rmax<10)
columns=c("regulator",level)
if(length(rid)>0){
topRegs1=topRegs[-rid,columns]
}else{
topRegs1=topRegs;
}


topRegs1[is.na(topRegs1)]=0
topRegs1$id=rowSums(topRegs1[,-1]>0)
topRegs2=topRegs1[order(topRegs1$id,decreasing = T),columns]
dtRegs1=topRegs2 %>% pivot_longer(2:ncol(topRegs2),names_to="cluster",values_to="ntarget")
dtRegs1$regulator=factor(dtRegs1$regulator,levels=topRegs2$regulator)
fwrite(dtRegs1,file=paste0(outpath,"/kmeans_top",rtop,"regulators_bubbleplot.txt"),sep="\t",row.names = F,quote=F,col.names = T)
write.table(topRegs2$regulator,file=paste0(outpath,"/kmeans_top",rtop,"regulators_level.txt"),sep="\n",row.names = F,quote=F,col.names = F)

#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
#library(ggpubr)
#library(scales)
my_cols <-c("#001B87FF","#31DB92FF","#FCC201FF")        # yellow-green-blue
dtRegs1=fread(paste0(outpath,"/kmeans_top5regulators_bubbleplot.txt"))
reglevel=topRegs2$regulator #fread("kmeans_top5regulators_level.txt",header = F)
outname=paste0(outpath,"/kmeans_top5regulators_bubbleplot.pdf")
#dtRegs1$cluster <- factor(dtRegs1$cluster,levels=columns)
#dtRegs1$regulator <- factor(dtRegs1$regulator,levels=rev(reglevel$V1))

dtRegs1$xaxis<- factor(dtRegs1$cluster,levels=columns)
#dtRegs1$yaxis<- factor(dtRegs1$regulator,levels=rev(reglevel$V1))
dtRegs1$yaxis<- factor(dtRegs1$regulator,levels=rev(reglevel))
pdf(outname, height = 10,width = 12,useDingbats=FALSE)
#ggballoonplot(dtRegs1, x="cluster", y="regulator", size="ntarget", fill = "ntarget")+
#  scale_fill_gradientn(colors = my_cols,limits=c(0,max(dtRegs1$ntarget)),oob=squish)
  cdat_sp <- ggplot(dtRegs1, aes(x = xaxis, y = yaxis, size = ntarget,color=ntarget)) +
    geom_point()+ 
    ggtitle(label = "top regulators")+
    theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(dtRegs1$ntarget)))+
   # theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
	theme(legend.position="bottom")+
    theme(text = element_text(size=11),axis.text = element_text(size = 11),axis.text.x = element_text(size = 11.5),
          strip.text = element_text(size = 11.5))#+ facet_wrap(~ncell,ncol=3, labeller = "label_both")

plot(cdat_sp);
dev.off()
#  rremove("grid")       # remove grid lines




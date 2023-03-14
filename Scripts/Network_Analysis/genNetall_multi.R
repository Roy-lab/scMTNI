genNetall_multi=function(prefix,data1,outname,w,h,nr=2,ntop=20,cellorder){
  # https://journal.r-project.org/archive/2017/RJ-2017-023/RJ-2017-023.pdf
  # https://cran.r-project.org/web/packages/ggnetwork/vignettes/ggnetwork.html#geom_edges
  library(dplyr)
  library(GGally)
  library(network)
  library(sna)
  library(ggplot2)
  library(tidyverse)
  library(gridExtra)
  library(ggnetwork)

  names(data1)=c("source","destination","weight","cell") ## edges
  data1$source=as.character(data1$source)
  data1$destination=as.character(data1$destination)
  ## get nodes:
  nodes=unique(c(data1$source,data1$destination))
  dnodes=data.frame(id=1:length(nodes),label=nodes) #nodesize=1
  edges <- data1 %>% left_join(dnodes, by = c("source" = "label")) %>% rename(from = id)
  edges <- edges %>% left_join(dnodes, by = c("destination" = "label")) %>% rename(to = id)
  edges <- select(edges, from, to, weight,cell)
  print(nrow(edges))
  dedges=unique(edges)
  dnodes$label=gsub("\\.","-",dnodes$label)
  mynetwork <- network(dedges, vertex.attr = dnodes, matrix.type = "edgelist", ignore.eval = FALSE,directed = T,multiple=T)
  # assign edge attributes (day)
  set.edge.attribute(mynetwork, "edgePresense", ifelse(mynetwork %e% "weight" > 0, "edge present", "edge absent"))
  #set.edge.attribute(mynetwork, "edgealpha", ifelse(mynetwork %e% "weight" > 0, 1, 0))

  ## SP Mar/14/2023: This line makes unexpected/gibbrish changes to the "mynetwork" object in
  ## version 1.18.1 of the "network" package. On the contrary, this line does
  ## not have any effect on the said object in version 1.17.1 of the "network"
  ## package. Either way, this line is not required. Hence, I'm commenting it.
  # set.edge.attribute(mynetwork, "cell", dedges[, 4])

  #set.edge.attribute(mynetwork, "edgename", dedges[, 5])
  set.seed(10052016)
  gnet=ggnetwork(mynetwork, arrow.size = 0.1,arrow.gap = 0.015, by = "cell",weights = "weight",
                 #layout = "segeo") #
                 layout = "fruchtermanreingold") #"fruchtermanreingold"
                 #layout = "kamadakawai") 

  ## get out-degree for each regulator node per cell:
  #cells=unique(data1$cell)
  cells=cellorder;
  dnode=NULL
  for(cell in cells){
    id=which(data1$cell==cell)
    dcell=data1[id,]
    dreg=dcell%>%group_by(source) %>% summarise(ntarget=n()) %>% arrange(desc(ntarget)) %>% as.data.frame()
    dreg$labelshort=as.character(dreg$source)
    if(nrow(dreg)>ntop){
      cf=dreg$ntarget[ntop]
      dreg$labelshort[which(dreg$ntarget<cf)]=""
    }##Add this to handle 0 edges for this node
    else if(nrow(dreg)==0){
	next
    }
    dreg$cell=cell
    dnode=rbind(dnode,dreg)
  }
  names(dnode)[1]="label"
  gnet=merge(gnet,dnode,by=c("label","cell"),all.x=T)
  gnet$ntarget[which(is.na(gnet$ntarget))]=0
  gnet$labelshort[which(is.na(gnet$labelshort))]=""  
  gnet$cell=gsub("cluster","C",gnet$cell)
  gnet$cell=factor(gnet$cell,levels = gsub("cluster","C",cells))
  cols <- c("edge present" = "salmon", "edge absent" = "grey")

  ## plot networks:
  #pdf(outname,width = 7,height = 6,useDingbats=FALSE)#pdf("test_MEF.pdf",width = 15,height = 10)
  #png(outname,width = w,height = h,res=300)
  #png(outname,width = w,height = h)
  #pdf(outname,width = w,height = h,useDingbats=FALSE)
  #pdf(outname,width = 20,height = 10,useDingbats=FALSE)
  #pdf(outname,width = 11,height = 8,useDingbats=FALSE)
  #pdf(outname,width = 11,height = 8,useDingbats=FALSE)
  pdf(outname,width = 11,height = 8,useDingbats=FALSE)
  #pdf(outname,useDingbats=FALSE)
  p=ggplot(gnet, aes(x, y, xend = xend, yend = yend)) +
    geom_edges(
      aes(color = edgePresense,alpha=weight),size=0.5, #size=weight+0.1
      arrow = arrow(length = unit(3, "pt"), type = "closed")) +
    geom_nodes(aes(size=ntarget),color = "lightskyblue",alpha=0.9) +  #aes(color = "lightskyblue"), 
    geom_nodetext(aes(label = labelshort,size=ntarget-1),color = "black")+  #fontface = "bold"
    theme(aspect.ratio=1,legend.position = "bottom")+
    scale_color_manual(values=cols)+ggtitle(prefix)+theme_void()+
  facet_wrap(~cell, nrow = nr) 
  #theme_facet(legend.position = "bottom",strip.text.x = element_text(size = 30,face="bold"))
  plot(p)
  dev.off()
}


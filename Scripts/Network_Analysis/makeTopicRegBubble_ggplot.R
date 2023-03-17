library(ggplot2)
library(data.table)

codename_prefix="Results_subsample/analysis/lda_TFcellbygene/k10/"

datatab1=NULL

## For regulator nodes with out-degree > "maxdeg,"
## reduce their out-degree to "maxdeg."
## It is done to upper bound the size of the
## bubbles corr. to these nodes so that
## the bubble plot looks more aesthetic.
maxdeg=40

outpath=codename_prefix

## For each topic
for(i in 1:10){
  datafile <- paste(c(codename_prefix,"topic",i,"_giantcomponent_min10_regdegree.txt"),collapse="")

  xordfile <- "ExampleData/celltype_order.txt"

  yordfile <- paste(c(codename_prefix,"topic",i,"_min10_yorder.txt"),collapse="")

  data <- read.csv(datafile,header=F,sep='\t')
  # datatab <- data.frame("yaxis"=data$V1, "xaxis"=data$V2, "Degree"=data$V3)
  datatab <- data.frame("yaxis"=data$V1, "xaxis"=data$V2, "ntarget"=data$V3)

  ## For regulator nodes with out-degree > "maxdeg,"
  ## reduce their out-degree to "maxdeg."
  ## It is done to upper bound the size of the 
  ## bubbles corr. to these nodes so that
  ## the bubble plot looks more aesthetic.
  # ids=which(datatab$Degree>maxdeg)
  ids=which(datatab$ntarget>maxdeg)
  # datatab$Degree[ids]=maxdeg;
  datatab$ntarget[ids]=maxdeg;

  xord <- read.csv(xordfile,header=F)
  datatab$xaxis <- factor(datatab$xaxis,levels=xord$V1)

  yord <- read.csv(yordfile,header=F)
  datatab$yaxis <- factor(datatab$yaxis,levels=rev(yord$V1))

  datatab$topic=i
  datatab1=rbind(datatab1,datatab)
}


##############################################################
## Begin: Generate a bubble plot suitable for the manuscript
###############################################################
# outname=paste0(outpath,"alltopics_rtop",rtop,"maxdeg",maxdeg,"_regulators_min10_bubbleplot_v2.pdf")
outname=paste0(outpath,"alltopics_","maxdeg",maxdeg,"_regulators_min10_bubbleplot.pdf")
pdf(outname, height = 6,width = 15,useDingbats=FALSE)

## yellow-green-blue
my_cols <-c("#001B87FF","#31DB92FF","#FCC201FF")

# cdat_sp <- ggplot(datatab1, aes(x = xaxis, y = yaxis, size = Degree,color=Degree)) +
cdat_sp <- ggplot(datatab1, aes(x = xaxis, y = yaxis, size = ntarget,color = ntarget)) +
  geom_point()+
  # ggtitle(label = "top regulators")+
  ggtitle(label = "Regulators with in+out degree >= 10")+
  xlab(element_blank())+
  ylab(element_blank())+
  # theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$Degree)))+
  theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$ntarget)))+
  # theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
  theme(legend.position="bottom")+
  # theme(text = element_text(size=8),axis.text = element_text(size = 8),axis.text.x = element_text(size = 6.5),
  theme(text = element_text(size=8),axis.text.y = element_text(size = 6.5),axis.text.x = element_text(size = 6.5, angle = 90),
        strip.text = element_text(size = 8))+scale_size(range = c(0.1,5))+facet_wrap(~topic,ncol=5, scales="free_y",labeller = "label_both")

plot(cdat_sp);
dev.off()
############################################################
## End: Generate a bubble plot suitable for the manuscript
###########################################################


#####################################################
## Begin: Generate a bubble plot suitable for slides
#####################################################
# outname=paste0(outpath,"alltopics_rtop",rtop,"maxdeg",maxdeg,"_regulators_min10_bubbleplot_v2.pdf")
outname=paste0(outpath,"alltopics_","maxdeg",maxdeg,"_regulators_min10_bubbleplot_v2.pdf")
# pdf(outname, height = 6,width = 15,useDingbats=FALSE)
pdf(outname, height = 8,width = 15,useDingbats=FALSE)

## yellow-green-blue
my_cols <-c("#001B87FF","#31DB92FF","#FCC201FF")

# cdat_sp <- ggplot(datatab1, aes(x = xaxis, y = yaxis, size = Degree,color=Degree)) +
cdat_sp <- ggplot(datatab1, aes(x = xaxis, y = yaxis, size = ntarget,color = ntarget)) +
  geom_point()+
  # ggtitle(label = "top regulators")+
  ggtitle(label = "Regulators with in+out degree >= 10")+
  xlab(element_blank())+
  ylab(element_blank())+
  # theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$Degree)))+
  theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$ntarget)))+
  # theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
  # theme(legend.position="bottom")+
  theme(legend.position="none")+
  # theme(text = element_text(size=8),axis.text = element_text(size = 8),axis.text.x = element_text(size = 6.5),
  # theme(text = element_text(size=8),axis.text.y = element_text(size = 6.5),axis.text.x = element_text(size = 6.5, angle = 90),
  theme(text = element_text(size=8),axis.text.y = element_text(size = 6.5),axis.text.x = element_text(size = 6.5, angle = 15),
        # strip.text = element_text(size = 8))+scale_size(range = c(0.1,5))+facet_wrap(~topic,ncol=5, scales="free_y",labeller = "label_both")
        strip.text = element_text(size = 8))+scale_size(range = c(0.1,5))+facet_wrap(~topic,ncol=4, scales="free_y",labeller = "label_both")

plot(cdat_sp);
dev.off()
#####################################################
## End: Generate a bubble plot suitable for slides
#####################################################



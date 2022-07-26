
#!/user/bin/env Rscript

library(rliger)
library(data.table)
fnames=c("scATACseq_FBS.txt",
         "scRNAseq_FBS.txt")

mypanel4=list()
samples=list('atac','rna');
for (i in 1:length(samples))
{
  print(paste0("Reading sample",i))
  fname=fnames[i]
  myCount<-fread(fname,header=T)
  #genenames<-read.table(paste(data_dir,"/common_genes.txt",sep=""))
  celllabel=colnames(myCount)[2:ncol(myCount)]
  if(i==2){
    barcode=do.call(rbind,strsplit(celllabel,"-"))
    cellbarcode=paste0(barcode[,1],"_",barcode[,2],"-",barcode[,3],"-",barcode[,4])
  }else{
    cellbarcode=celllabel
  }
  
  myCount=as.data.frame(myCount)
  genenames=myCount[,1]
  rownames(myCount)<-genenames
  print("Added gene names")
  #barcodes<-read.table(paste(data_dir,"/allbarcodes_",samples[i],".txt",sep=""))
  #colnames(myCount)<-barcodes[,1]
  print("Added colnames")
  suff=samples[[i]]
  myCount=myCount[,-1]
  colnames(myCount)<-paste(suff,cellbarcode,sep="_")
  #mypanel4[[suff]]=as(myCount,"dgCMatrix")
  mypanel4[[suff]]=myCount
}

# atac: 12465 53883
# rna: 12465  5172 cells
print("Done reading all samples")
outpath="/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/liger/sqrt_atacsum100_rnacell50_FBS/"
dir.create(outpath)
setwd(outpath)
save(mypanel4,genenames,file="liger_data.Rdata")

## Stage I: Preprocessing and Normalization 
ligerdata=createLiger(mypanel4)
ligerdata=normalize(ligerdata)
#Note that by setting the parameter datasets.use to 2, genes will be selected only from the scRNA-seq dataset (the second dataset) by the selectGenes function. 
##We recommend not using the ATAC-seq data for variable gene selection because the statistical properties of the ATAC-seq data are very different from scRNA-seq, 
###violating the assumptions made by the statistical model we developed for selecting genes from RNA data. 
#varthre=0.05
pdf(paste0("variablegenes.pdf"))
ligerdata=selectGenes(ligerdata,datasets.use = 2,do.plot=T) # var.thresh = varthre (0.1 default)
#sr_liger = selectGenes(sr_liger, var.thresh = 1)
#ligerdata@var.genes=genenames;
print(length(ligerdata@var.genes))
dev.off()
##Finally, the scaleNotCenter function scales normalized datasets without centering by the mean, giving the nonnegative input data required by iNMF.
ligerdata = scaleNotCenter(ligerdata)
norm_rna=as.data.frame(as.matrix(ligerdata@norm.data[["rna"]]))
norm_atac=as.data.frame(as.matrix(ligerdata@norm.data[["atac"]]))
fwrite(norm_rna,file=paste0("scRNAseq_liger_normliazed.txt"),quote=F,sep="\t",row.names = T,col.names = T)
fwrite(norm_atac,file=paste0("scATACseq_liger_normliazed.txt"),quote=F,sep="\t",row.names = T,col.names = T)


plotByDatasetAndCluster <- function(object, clusters = NULL, title = NULL, pt.size = 0.3,
                                    text.size = 3, do.shuffle = TRUE, rand.seed = 1,
                                    axis.labels = NULL, do.legend = TRUE, legend.size = 5,
                                    reorder.idents = FALSE, new.order = NULL,
                                    return.plots = FALSE, legend.fonts.size = 12) {
  library(dplyr)
  tsne_df <- data.frame(object@tsne.coords)
  colnames(tsne_df) <- c("tsne1", "tsne2")
  tsne_df[['Dataset']] <- unlist(lapply(1:length(object@H), function(x) {
    rep(names(object@H)[x], nrow(object@H[[x]]))
  }))
  if (reorder.idents == TRUE){
    tsne_df$Dataset <- factor(tsne_df$Dataset, levels = new.order)
  }
  c_names <- names(object@clusters)
  if (is.null(clusters)) {
    # if clusters have not been set yet
    if (length(object@clusters) == 0) {
      clusters <- rep(1, nrow(object@tsne.coords))
      names(clusters) <- c_names <- rownames(object@tsne.coords)
    } else {
      clusters <- object@clusters
      c_names <- names(object@clusters)
    }
  }
  tsne_df[['Cluster']] <- clusters[c_names]
  
  celltype=do.call(rbind,strsplit(c_names,"-"))
  tsne_df[['Celltype']] <- celltype[,1]#,paste0(celltype[,1],"-",celltype[,2])
  if (do.shuffle) {
    set.seed(rand.seed)
    idx <- sample(1:nrow(tsne_df))
    tsne_df <- tsne_df[idx, ]
  }
  
  p1 <- ggplot(tsne_df, aes_string(x = 'tsne1', y = 'tsne2', color = 'Dataset')) + theme_bw() +
    theme_cowplot(legend.fonts.size) + geom_point(size = pt.size, stroke = 0.2) +
    guides(color = guide_legend(override.aes = list(size = legend.size)))
  
  centers <- tsne_df %>% group_by(.data[['Cluster']]) %>% summarize(
    tsne1 = median(x = .data[['tsne1']]),
    tsne2 = median(x = .data[['tsne2']]),
  )
  p2 <- ggplot(tsne_df, aes_string(x = 'tsne1', y = 'tsne2', color = 'Cluster')) +
    theme_cowplot(legend.fonts.size) + geom_point(size = pt.size, stroke = 0.2) +
    geom_text(data = centers, mapping = aes_string(label = 'Cluster'), colour = "black", size = text.size) +
    guides(color = guide_legend(override.aes = list(size = legend.size)))
  
  centers1 <- tsne_df %>% group_by(.data[['Celltype']]) %>% summarize(
    tsne1 = median(x = .data[['tsne1']]),
    tsne2 = median(x = .data[['tsne2']]),
  )
  
  p3 <- ggplot(tsne_df, aes_string(x = 'tsne1', y = 'tsne2', color = 'Celltype')) +
    theme_cowplot(legend.fonts.size) + geom_point(size = pt.size, stroke = 0.2) +
    geom_text(data = centers1, mapping = aes_string(label = 'Celltype'), colour = "black", size = text.size) +
    guides(color = guide_legend(override.aes = list(size = legend.size)))
  
  if (!is.null(title)) {
    p1 <- p1 + ggtitle(title[1])
    p2 <- p2 + ggtitle(title[2])
    p3 <- p3 + ggtitle(title[3])
  }
  if (!is.null(axis.labels)) {
    p1 <- p1 + xlab(axis.labels[1]) + ylab(axis.labels[2])
    p2 <- p2 + xlab(axis.labels[1]) + ylab(axis.labels[2])
    p3 <- p3 + xlab(axis.labels[1]) + ylab(axis.labels[2])
  }
  p1 <- p1 + theme_cowplot(12)
  p2 <- p2 + theme_cowplot(12)
  p3 <- p3 + theme_cowplot(12)
  if (!do.legend) {
    p1 <- p1 + theme(legend.position = "none")
    p2 <- p2 + theme(legend.position = "none")
    p3 <- p3 + theme(legend.position = "none")
  }
  if (return.plots) {
    return(list(p1, p2,p3))
  } else {
    print(p1)
    print(p2)
    print(p3)
  }
}
library(gridExtra)
library(ggplot2)
for(k in c(8,10,12,15,20)){  
  outdir=paste0(outpath,"/k",k,"/")
  dir.create(outdir,recursive = T)
  load(paste0(outdir,"liger_k",k,".Rdata"))
  # Stage II: Joint Matrix Factorization
  # perform joint matrix factorization (iNMF) on the normalized and scaled RNA and ATAC data. 
  # This step calculates metagenes–sets of co-expressed genes that distinguish cell populations–containing both shared and dataset-specific signals. 
  # The cells are then represented in terms of the “expression level” of each metagene, providing a low-dimensional representation that can be used 
  # for joint clustering and visualization. To run iNMF on the scaled datasets, we use the optimizeALS function with proper hyperparameter settings.
  ## (1) k: We find that a value of k in the range 20 - 40 works well for most datasets.
  ## (2) lambda: a regularization parameter. Larger values penalize dataset-specific effects more strongly, causing the datasets to be better aligned,
  ## but possibly at the cost of higher reconstruction error. The default value is 5.
  ## (3) thresh: This sets the convergence threshold. Lower values cause the algorithm to run longer. The default is 1e-6.
  ## (4) max.iters: This variable sets the maximum number of iterations to perform. The default value is 30.
  sr_liger=optimizeALS(ligerdata,k=k,print.obj=T)
  # tsne_before<- runTSNE(sr_liger, use.raw = T)
  # pdf(paste0(outdir,"/TSNE_unalignedcellfactor.pdf"),height=5,width=10,useDingbats=FALSE)
  # p=plotByDatasetAndCluster(tsne_before,return.plots = T)
  # grid.arrange(p[[1]],p[[2]], ncol=2)
  # dev.off()
  
  # Stage III: Quantile Normalization and Joint Clustering (1 minute)
  # Using the metagene factors calculated by iNMF, we assign each cell to the factor on which it has the highest loading, giving joint clusters that correspond across datasets. 
  # We then perform quantile normalization by dataset, factor, and cluster to fully integrate the datasets.
  sr_liger=quantile_norm(sr_liger)
  #sr_liger <- louvainCluster(sr_liger, resolution = 0.2)
  tsne_after=runTSNE(sr_liger)
  #pdf(paste0(outdir,"/TSNE_jointclustered.pdf"),height=5,width=15,useDingbats=FALSE)
  #png(paste0(outdir,"/TSNE_jointclustered.png"),height=480*3,width=480*12,res=300)
  pdf(paste0(outdir,"/TSNE_jointclustered1.pdf"),height=5,width=10,useDingbats=FALSE)
  p=plotByDatasetAndCluster(tsne_after,return.plots = T)
  grid.arrange(p[[2]],p[[3]],ncol=2)
  #grid.arrange(p[[1]],p[[2]],p[[3]], layout_matrix = matrix(c(1,2,3,3),1),ncol=3)
  dev.off()
  
  pdf(paste0(outdir,"/plotofclusters_sizedby_proportionof_totalcells.pdf"),useDingbats=FALSE)
  plotClusterProportions(sr_liger)
  dev.off();
  
  pdf(paste0(outdir,"/scatterplots_unaligned_and_aligned_factor_loadings.pdf"),useDingbats=FALSE)
  plotFactors(sr_liger)
  dev.off();
  #  Plots heatmap of matrix, with red representing high total loadings for a factor, black low
  pdf(paste0(outdir,"/heatmap_cluster_correspondence.pdf"),useDingbats=FALSE)
  plotClusterFactors(sr_liger,use.aligned=T)
  dev.off();
  
  
  pumap <- runUMAP(sr_liger,n_neighbors = 30) #distance = 'cosine',  min_dist = 0.3
  #png(paste0(outdir,"/UMAP_jointclustered_euclidean_n30.png"),height=480,width=480*3.5,useDingbats=FALSE)
  png(paste0(outdir,"/UMAP_jointclustered_euclidean_n30.png"),height=480*3,width=480*12,res=300)
  p=plotByDatasetAndCluster(pumap, axis.labels = c('UMAP 1', 'UMAP 2'),return.plots = T)
  grid.arrange(p[[1]],p[[2]],p[[3]], layout_matrix = matrix(c(1,2,3,3),1),ncol=3)
  dev.off();
  
  pumap <- runUMAP(sr_liger,n_neighbors = nb) #distance = 'cosine',  min_dist = 0.3
  pdf(paste0(outdir,"/UMAP_jointclustered_euclidean_n",nb,".pdf"),height=5,width=10,useDingbats=FALSE)
  p=plotByDatasetAndCluster(pumap, axis.labels = c('UMAP 1', 'UMAP 2'),return.plots = T)
  grid.arrange(p[[2]],p[[3]],ncol=2) #layout_matrix = matrix(c(1,2,2),1),
  dev.off();
  
  
  cfile=paste0(outdir,"/ligerclusters.txt")
  write.table(sr_liger@clusters,cfile,sep="\t",quote=F,col.names=FALSE)
  
  cfile=paste0(outdir,"/W.txt")
  write.table(sr_liger@W,cfile,sep="\t",quote=F,row.names=FALSE)
  cfile=paste0(outdir,"/Hnorm.txt")
  write.table(sr_liger@H.norm,cfile,sep="\t",quote=F,col.names=FALSE)
  save(sr_liger,file=paste0(outdir,"liger_k",k,".Rdata"))
  
  ## April 15, 2021: cell distribution:
  #outdir=paste0(outpath,"/k",k,"/")
  #dir.create(outdir,recursive = T)
  #load(paste0(outdir,"liger_k",k,".Rdata"))
  plotdistribution(sr_liger,k)
  plotRNAdistribution(sr_liger,k)
}

plotdistribution=function(sr_liger,k){
  prefix=paste0("liger_k",k)
  barcode=do.call(rbind,strsplit(names(sr_liger@clusters),"-"))
  cluster0=sr_liger@clusters
  cellbarcode=barcode[,1]#,paste0(barcode[,1]"-",barcode[,2])
  cluster=data.frame(cluster=cluster0,cellbarcode)
  dt=table(cluster[c("cluster","cellbarcode")])
  #dt1=data.frame(dt/rowSums(dt))
  totalcell=rowSums(dt)
  
  dt0=data.frame(cluster=names(totalcell),cellbarcode="Totalcells",Freq=totalcell,Totalcells=1)
  dt=as.matrix(dt)
  atac=dt[,grep("atac",colnames(dt))]
  rna=dt[,grep("rna",colnames(dt))]
  totalcell=rowSums(atac)
  dt0=rbind(dt0,data.frame(cluster=names(totalcell),cellbarcode="atac",Freq=totalcell,Totalcells=1))
  atac1=data.frame(atac/rowSums(atac))
  rna1=data.frame(rna/rowSums(rna))
  dt1=rbind(atac1,rna1)
  dt1$Totalcells=0
  totalcell=rowSums(rna)
  dt0=rbind(dt0,data.frame(cluster=names(totalcell),cellbarcode="rna",Freq=totalcell,Totalcells=1))
  dt1$cluster=paste0("C",dt1$cluster)
  dt0$cluster=paste0("C",dt0$cluster)
  mylevels=paste0("C",1:k)
  
  library(ggplot2)
  setwd(outdir)
  #outname=paste0(prefix,"_celldistribution.pdf")
  
  #pdf(outname, height = 4,width =20)
  p1=ggplot(dt0, aes(cellbarcode, cluster)) + geom_tile(aes(fill = Freq),colour = "white") + geom_text(aes(label = round(Freq,digits = 2)))+
    scale_x_discrete(position = "top")+ylim(rev(mylevels))+ggtitle(paste0(prefix," clusters"))+
    scale_fill_gradient(low = "white",high = "blue")+ theme_bw()+
    theme(text = element_text(size=8),legend.position="bottom")
  p2=ggplot(dt1, aes(cellbarcode, cluster)) + geom_tile(aes(fill = Freq),colour = "white") + geom_text(aes(label = round(Freq,digits = 2)))+
    scale_x_discrete(position = "top")+ylim(rev(mylevels))+ggtitle(paste0(prefix," clusters"))+
    scale_fill_gradient(low = "white",high = "red")+ theme_bw()+
    theme(text = element_text(size=8),legend.position="bottom") #axis.text.x  = element_text(angle=90, vjust=0.5, size=8)
  
  library(gridExtra)
  # grid.arrange(p1,p2, ncol=2, layout_matrix = matrix(c(1,2,2,2),nrow=1) )
  # dev.off()
  
  dtc=table(cluster[c("cellbarcode","cluster")])
  dtc3=data.frame(dtc/rowSums(dtc))
  totalcell=rowSums(dtc)
  dtc3$Totalcells=0
  dtc2=data.frame(samples=names(totalcell),cellbarcode="Totalcells",Freq=totalcell,Totalcells=1)
  dtc=as.matrix(dtc)
  dtc3$cluster=paste0("C",dtc3$cluster)
  dtc3$cluster=factor(dtc3$cluster,levels = paste0("C",1:k))
  #dtc2$cluster=paste0("cluster",dtc2$cluster)
  mylevels=levels(dtc3$cellbarcode)
  outname=paste0(prefix,"_alldistributions.pdf")
  
  
  p3=ggplot(dtc2, aes(cellbarcode, samples)) + geom_tile(aes(fill = Freq),colour = "white") + geom_text(aes(label = round(Freq,digits = 2)))+
    scale_x_discrete(position = "top")+ylim(rev(mylevels))+ggtitle(paste0(prefix," clusters"))+
    scale_fill_gradient(low = "white",high = "blue")+ theme_bw()+
    theme(text = element_text(size=8),legend.position="bottom")
  p4=ggplot(dtc3, aes(cluster,cellbarcode)) + geom_tile(aes(fill = Freq),colour = "white") + geom_text(aes(label = round(Freq,digits = 2)))+
    scale_x_discrete(position = "top")+ylim(rev(mylevels))+ggtitle(paste0(prefix," clusters"))+
    scale_fill_gradient(low = "white",high = "red")+ theme_bw()+
    theme(text = element_text(size=8),legend.position="bottom")
  
  library(gridExtra)
  pdf(outname, height = 4,width =12,useDingbats=FALSE)
  grid.arrange(p1,p2,p3,p4, ncol=4, layout_matrix = matrix(c(1,2,2,3,4,4),nrow=1) )
  dev.off()
}


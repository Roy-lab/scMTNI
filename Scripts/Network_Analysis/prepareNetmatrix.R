#!/user/bin/env Rscript

## Manage libraries -----
library(data.table)
library(Matrix)
library(stringr)

## Manage input args -----
args=commandArgs(trailingOnly = T);
if(length(args)!= 5 && argslength(args) != 6)
{
  stop("Exactly 5 or 6 arguments must be supplied (cf, n, cellfile, indir, outdir, [binary]).\n", call. = FALSE)
}else
{
  print("Start!")
}

#### OLD ARGUMENTS
    #cellfile=args[1]  ## celltype_order.txt
    #indir=args[2]  ## Results_/subsample/analysis/
    #filter=args[3] #"_full" or "_cf0.8"
####

cf=args[1] #"top" or "cf"
n=as.numeric(args[2]) # eg n=5000 for top edge of cf=0.8 for cf  
cellfile=args[3] # celltype_order.txt
indir=args[4] # input path
outdir=args[5]  # output path 

if(length(args)==5){
  binary=""
}else{
  binary=args[6]
}


### Error checking -----
#### error checking arg1 -----
if (cf != "top" && cf != "cf") {
  stop("Invalid value for cf. Use 'top' or 'cf'.", call. = FALSE)
}

#### error checking arg2 -----
if (cf == "top")
{
  if(n <= 0) {
    stop("n must be a positive integer when cf is 'top'.", call. = FALSE)
  }
  print(paste("Extracting top", n, "edges for each cell type."))
} else if (cf == "cf") {
  if(n < 0 || n > 1) {
    stop("n must be a number between 0 and 1 when cf is 'cf'.", call. = FALSE)
  }
  print(paste("Extracting edges with confidence threshold of", n, "for each cell type."))
}


#### error checking arg3 -----
if (file.exists(cellfile) == FALSE) {
  stop(paste("Cell file", cellfile, "does not exist."), call. = FALSE)
}


#### error checking arg4 -----
if(dir.exists(indir) == FALSE) {
  stop(paste("Input path", indir, "does not exist."), call. = FALSE)
}

#### error checking arg5 -----
if(dir.exists(outdir) == FALSE) {
  print(paste("Output path", outdir, "does not exist."))
  dir.create(outdir, recursive = TRUE)
  print(paste("Created output path:", outdir))
}

#### error checking arg6 -----
if(binary == "")
{
  print("Binary option not provided. Defaulting to no binary conversion.")
}

### Read order file -----
cells=read.table(cellfile,stringsAsFactors = F)
cells=cells$V1

genes=NULL
regulators=NULL
for(cell in cells){
  
  ### Read cell data file -----
  infile <- paste0(indir,'/', cell, "/consensus_edges_", cf, n, ".txt")
  if(file.exists(infile) == FALSE) {
    stop(paste("Input file for cell type", cell, "does not exist:", infile), call. = FALSE)
  }
  d=fread(infile)
  
  ### Apply filter ------
  d=d[order(d$V3,decreasing = T),]
  
  #### Apply topk filter -----
  if(cf=="top"){
    #d1=d[1:n,] THIS IS PROBLEMATIC SINCE IT DOES NOT HANDLE TIES! SHS
    
    # Get the n-th largest value, or last value if n exceeds network size 
    threshold <- sort(d$V3, decreasing = TRUE)[min(n, nrow(d))] 
    
    # apply top n_filter,
    id <- which(d$V3 >= threshold)  # Select rows with value >= n-th largest
    
  #### Apply cf filter -----
  }else if(cf=="cf"){
    id <- which(d$V3>=n)
  
  }else{
    stop("Invalid value for cf. Use 'top' or 'cf'.", call. = FALSE)
  }
  d = d[id, ]
  
  ### Identify unique regulators -----
  regulators=c(regulators,unique(d$V1))
  reg_count=length(regulators)
  print(paste0("The number of regulators for ",cell," are ",reg_count))
  genes=c(genes,unique(d$V2))
}

### Identify unique genes -----
genes=unique(genes)

### Generate genes and regulator data frames -----
genes=data.frame(gene=genes,id=1:length(genes))
genes$gene=as.character(genes$gene)
write.table(genes,file=paste0(outdir,"genes_filteredlowexpression",filter,binary,".txt"),quote=F,row.names=F,col.names=F,sep="\t")


regulators=unique(regulators)
regulators=data.frame(gene=regulators)
regulators=merge(regulators,genes,by="gene")
regulators=regulators[order(regulators$id),]
regid=regulators$id
write.table(regid,file=paste0(outdir,"TFtarget_regulator_id",filter,binary,".txt"),quote=F,row.names=F,col.names=F,sep="\n")
regulators$id=1:nrow(regulators)
regulators$gene=as.character(regulators$gene)


for(cell in cells){
  infile <- paste0(indir,'/', cell, "/consensus_edges_", cf, n, ".txt")
  if(file.exists(infile) == FALSE) {
    stop(paste("Input file for cell type", cell, "does not exist:", infile), call. = FALSE)
  }
  d=fread(infile)
  
  ### Apply filter ------
  d=d[order(d$V3,decreasing = T),]
  
  #### Apply topk filter -----
  if(cf=="top"){
    #d1=d[1:n,] THIS IS PROBLEMATIC SINCE IT DOES NOT HANDLE TIES! SHS
    
    # Get the n-th largest value, or last value if n exceeds network size 
    threshold <- sort(d$V3, decreasing = TRUE)[min(n, nrow(d))] 
    
    # apply top n_filter,
    id <- which(d$V3 >= threshold)  # Select rows with value >= n-th largest
    
    #### Apply cf filter -----
  }else if(cf=="cf"){
    id <- which(d$V3>=n)
    
  }else{
    stop("Invalid value for cf. Use 'top' or 'cf'.", call. = FALSE)
  }
  d = d[id, ]
  
  regs=data.frame(gene=unique(d$V1))
  regs=merge(regulators,regs,by="gene")
  regs=regs[order(regs$id),]
  #print(paste0("The number of regulators for ",cell," are ",regs,"\n"))
  #regs$id=1:nrow(regs)
  #dmat=matrix(0,nrow(regs),nrow(genes))
  dmat=matrix(0,nrow(regulators),nrow(genes))
  for(i in regs$id) #1:nrow(regs))
  {
    #reg=regs$gene[i]
    reg=regulators$gene[i]
    d1=d[which(d$V1==reg),]
    d1=merge(d1,genes,by.x="V2",by.y="gene")
    dmat[i,d1$id]=d1$V3
  }
  #dmat=forceSymmetric(dmat, "UP")
  dmat=data.table(as.matrix(dmat))
  row.names(dmat)=regulators$gene  #regs$gene 
  names(dmat)=genes$gene
  fwrite(dmat,file=paste0(outdir,cell,"_consensus_edges_filteredlowexpression_mat",filter,binary,".txt"),quote=F,row.names=T,col.names=T,sep="\t")
}



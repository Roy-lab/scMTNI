#!/user/bin/env Rscript
##This is Shilu's script with SR updated to take an output path 
#install.packages("stringr", repos='http://cran.us.r-project.org')
library(stringr)

## Input args ----
args=commandArgs(trailingOnly = T);
if(length(args)==0)
{
  stop("At least one argument must be supplied(inputfile).n",call. = FALSE)
}else if(length(args)>0)
{
  print("Start!")
}

library(data.table)
cf=args[1] #"top" or "cf"
n=as.numeric(args[2]) # eg n=5000 for top edge of cf=0.8 for cf  
cellfile=args[3] # celltype_order.txt
inpath=args[4] # input path
outpath=args[5]  # output path 

## Error checking ----- 

### error checking arg1 ----- 
if (cf != "top" && cf != "cf") {
  stop("Invalid value for cf. Use 'top' or 'cf'.", call. = FALSE)
}


### error checking arg2 ----- 
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

### error checking arg3 ----- 
if (file.exists(cellfile) == FALSE) {
  stop(paste("Cell file", cellfile, "does not exist."), call. = FALSE)
}

### error checking arg4 -----
if(dir.exists(inpath) == FALSE) {
  stop(paste("Input path", inpath, "does not exist."), call. = FALSE)
}

### error checking arg5 -----
if(dir.exists(outpath) == FALSE) {
  print(paste("Output path", outpath, "does not exist."))
  dir.create(outpath, recursive = TRUE)
  print(paste("Created output path:", outpath))
}

## Read order file -----
cellsf=read.table(cellfile,stringsAsFactors = F)
cells=cellsf$V1

## Read data file -----

data=NULL
for(cell in cells){
  print(cell)
  
  ### Check if the input file exists -----
  infile=paste0(inpath,"/",cell,"/", "consensus_edges.txt")
  if (!file.exists(infile)) {
    stop(paste("Input file for cell type", cell, "does not exist:", infile), call. = FALSE)
  }
  d=fread(paste0(inpath,cell,"/consensus_edges.txt"))
  
  ### Apply filter ------
  d=d[order(d$V3,decreasing = T),]
  
  #### Apply topk filter -----
  if(cf=="top"){
    #d1=d[1:n,] THIS IS PROBLEMATIC SINCE IT DOES NOT HANDLE TIES! SHS
    
    # Get the n-th largest value, or last value if n exceeds network size 
    threshold = sort(d$V3, decreasing = TRUE)[min(n, nrow(d))] 
    
    # apply top n_filter,
    d1 <- d[d$V3 >= threshold, ]  # Select rows with value >= n-th largest
  
  #### Apply cf filter -----
  }else if(cf=="cf"){
    id=which(d$V3>=n)
    d1=d[id,]
  }
  else{
    stop("Invalid value for cf. Use 'top' or 'cf'.", call. = FALSE)
  }
  
  ### Write cell-type specific consensus edges -----
  print(paste(cell,nrow(d1)))
  fwrite(d1,file=paste0(outpath,"/",cell,"/consensus_edges_",cf,n,".txt"),quote=F,sep="\t",row.names = F,col.names = F)
  print(paste0(outpath,"/",cell,"/consensus_edges_",cf,n,".txt"))
  names(d1)=c("regulator","target",cell)
  
  ### Add cell-type specific consensus edges to the consensus data frame -----
  if(is.null(data)){
    data=d1
  }else{
    data=merge(data,d1,by=c("regulator","target"),all=T)
  }
}

#### Fill cell-type specific consensus edges to the consensus data frame -----
data[is.na(data)]=0
data$edge=paste0(data$regulator,"->",data$target)


### Write cell type-specific consensus edges to a single file -----
selcols=c("edge",cells)
dall=data[,..selcols]
fwrite(dall,file=paste0(outpath,"/consensus_edges_merged_",cf,n,".txt"),quote=F,sep="\t",row.names = F,col.names = T)



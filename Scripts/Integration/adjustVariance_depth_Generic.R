#' Adjust feature variance to remove effect of varying feature means.
#' Adapted from pagoda2: https://github.com/hms-dbmi/pagoda2
#'
#' @param counts Input data matrix
#' @param gam.k Number of additive models to use when fitting mean variance relationship
#' @param plot Whether or not to plot the mean variance relationship
#' @param max.adjusted.variance Maximum adjusted feature variance
#' @param min.adjusted.variance Minimum adjusted feature variance
#' @param verbose Whether or not to print output
#' @param q.val Maximum adjusted p-value for a feature to qualify as overdispersed
#'
#' @return Dataframe with adjusted feature variances
#'
#' @importFrom mgcv gam 
#' @import Matrix
#' @export
#'
AdjustVariance <- function(counts, gam.k = 10, plot = F, max.adjusted.variance = 1e3, min.adjusted.variance = 1e-3,
                           verbose = T, q.val = 0.05, gene_cnt_thresh=49, cell_cnt_thresh=2000,depth_norm=2000) {
  if (is.data.frame(counts) | is.matrix(counts)) {
    counts <- as(as.matrix(counts), "dgCMatrix")
  }
  counts <- Matrix::t(counts)
	
  cat("Depth normalizing with gene_th=",gene_cnt_thresh," cell_th=",cell_cnt_thresh, "\n")
  counts <- depthNormalize(counts,gene_cnt_thresh,cell_cnt_thresh,depth_norm)
  cat("Done depth normalizing\n")
  if(verbose) cat("calculating variance fit ...")
  df <- colMeanVarS(counts)
  #SR fix, Apr 6th 2020 to make sure counts is of the same dim as nz
  counts<-counts[,df$nz]
  cat("Returned from colMeanVarS\n") 		
  #df$oldm <- df$m; df$oldv <- df$v;
  df$m <- log(df$m); df$v <- log(df$v);
  #df$oldm <- log(df$oldm); df$oldv <- log(df$oldv);
  #cat(paste("colname length",length(colnames(counts)),"\n"))
  #cat(paste("colname length of df ",length(rownames(df)),"\n"))
  #rownames(df) <- colnames(counts);
  vi <- which(is.finite(df$v) & df$nobs >= 0);
  cat(length(vi))
  if(length(vi) < gam.k*1.5) { gam.k = 1 };# too few genes
  if(gam.k < 2) {
    if(verbose) cat(" using lm ")
    m <- lm(v ~ m, data = df[vi,])
  } else {
    if(verbose) cat(" using gam ")
    m <- mgcv::gam(as.formula(paste0("v ~ s(m, k = ", gam.k, ")")), data = df[vi,])
  }
  cat("done gamming")	
  df$res <- -Inf;  
  #getting the residual
  df$res[vi] <- resid(m,type='response')
  n.obs <- df$nobs; #diff(counts@p)
  df$lp <- as.numeric(pf(exp(df$res),n.obs,n.obs,lower.tail=F,log.p=F))
  df$lpa <- p.adjust(df$lp, method = "BH")
  n.cells <- nrow(counts)
  if(verbose) cat("calculating qv ... n.cells=",n.cells)
  df$qv <- as.numeric(qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=F)/n.cells)

  ods <- which(df$lpa < q.val)

  df$gsf <- geneScaleFactors <- sqrt(pmax(min.adjusted.variance,pmin(max.adjusted.variance,df$qv))/exp(df$v));
  df$gsf[!is.finite(df$gsf)] <- 0;
  odgenes <- rownames(df)[ods]
  df$overdispersed <- rownames(df) %in% odgenes
  allresults<-list("varinfo"=df,"normcnts"=counts)
  return(allresults)
}

depthNormalize <-function(counts,gene_cnt_thresh=49, cell_cnt_thresh=2000,depth_norm=2000){
##SR made updates Jan9, 2020. Using gene_cnt_thresh, cell_cnt_thresh and depth_norm as new params
	#thresh<-49
	thresh<-gene_cnt_thresh
	gsums<-Matrix::colSums(counts)
	#Here we are attempting to remove genes with <thresh cells expression. Genes are on columns
	cnames<-colnames(counts);
	cat("Got genenames names\n")
	counts<-counts[,gsums>thresh];
	nz<-gsums>thresh
	cat("select nz\n")
	colnames(counts)<-cnames[nz];
	#cat(paste("selecting",length(nz), "colnames\n"))
	cat(paste("selecting",length(which(nz)), "colnames\n"))

	#Added Sep 28th 2019, to remove cells with too few reads
	rnames<-rownames(counts)
	#thresh_c=5000;
	thresh_c=cell_cnt_thresh;
	depth <- Matrix::rowSums(counts)
	nz<-depth>thresh_c
	counts<-counts[nz,];
	depth<-depth[nz]
	rownames(counts)<-rnames[nz];
	counts<- counts/(depth/depth_norm)
	return(counts)
}

colMeanVarS <- function(counts){
	#get the sum
	print("in colMeanVar S")
	s_d<-colSums(counts)
	cat(length(s_d))
	nz<-s_d>0
	#SR change Apr 6th, 2020
	#s_d<-s_d[s_d>0]
	s_d<-s_d[nz]
	cat("Length after filtering ",length(s_d),"\n")
	#Calculate the mean
	cellcnt<-counts@Dim[1]
	print("Calculating mean")
	m1=(s_d/cellcnt)
	#Repeat the mean so we can turn this into  matrix
	v<-rep(m1,cellcnt)
	#genecnt<-counts@Dim[2]
	genecnt<-length(s_d)
	#reshape into  a matrix
	ms<-matrix(v,ncol=genecnt,byrow=TRUE)
	counts<-counts[,nz];
	diff<-(counts-ms)
	vsum<-diff*diff
	#get the variance
	print("Calculating variance")
	vars<-colSums(vsum)/cellcnt
	ncnts<-vector('integer',genecnt)
	#getting the number of obs
	print("Calculating nobs")
	#for (g in 1:genecnt)
	#{
	#	ncnts[g]<-sum(counts[,g]>0)
	#}
	ncnts<-colSums(counts>0)
	
	df<-data.frame("m"=m1,"v"=vars,"nobs"=ncnts,"nz"=which(nz))
  	rownames(df) <- colnames(counts);
	print('All done with meanvar calc')
	return(df)
}

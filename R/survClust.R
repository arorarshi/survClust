
##########################
# survClust
###########################

survclust<-function(combine.dist,survdat,k, cmd.k=NULL){
  if(is.null(rownames(survdat)))
    stop("rowanmes of survdat can't be NULL")
  
  if(is.null(rownames(combine.dist)))
    stop("rowanmes of combine.dist can't be NULL")
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("pdist", quietly = TRUE)) {
    stop("pdist package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  inter <- intersect(rownames(survdat), rownames(combine.dist))
  
  clin<-survdat[inter,]
  clin<-apply(clin,2,as.numeric)
  
  #run cmdscale
  combine.dist <- combine.dist[inter,inter]
  
  if(is.null(cmd.k)){cmd.k =nrow(combine.dist)-1 }
  if(!(is.null(cmd.k))){cmd.k =as.numeric(cmd.k) }
  
  cmd.combine.dist<-cmdscale(combine.dist,k=cmd.k, add=TRUE)$points
  
  my.k = as.numeric(k)
  survobj<-Surv(clin[,1], clin[,2])
  
  #run kmeans with 100 starts
  fit<-kmeans(cmd.combine.dist,my.k,nstart=100)
  
  #caluclate logrank
  fit.lr<-survdiff(Surv(clin[,1], clin[,2]) ~ fit$cluster)$chisq
  fit$fit.lr<-fit.lr
  
  #return fit and its logrank
  return(fit)
}

#predict test labels on survclust fitted
predict.test.label<-function(all.cmd,fit,k){
  all.cmd = as.matrix(all.cmd)
  train.snames = names(fit$cluster)
  test.snames = setdiff(rownames(all.cmd),train.snames)
  
  #where row - samples, col - genes
  centroid = matrix(NA, nrow = k, ncol = ncol(all.cmd))
  for (kk in 1:k) {
    #meaning k clust has one sample. #WARNING #check
    if(is.vector(all.cmd[names(fit$cluster)[which(fit$cluster==kk)],]) & ncol(all.cmd) > 1){
      message(paste0("k=",k, " training cluster has one sample, prediction might be inaccurate"))
      centroid[kk, ]=all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ]
    }
    
    if (!(is.null(dim(all.cmd[names(fit$cluster)[fit$cluster==kk],])))){
      if(ncol(all.cmd)> 1){centroid[kk, ]=apply(all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ], 2, mean)}
    }
    
    if(ncol(all.cmd)==1){centroid[kk,] = mean(all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ])}
  }
  
  dist.whole = apply(centroid,1,function(x) as.matrix(pdist(x,all.cmd)))
  
  #assign the cluster membership
  dist.labels = apply(dist.whole,1,which.min)
  names(dist.labels) = rownames(all.cmd)
  test.labels = dist.labels[test.snames]
  
  #is missing a class label via pdist
  if(length(unique(test.labels)) != k){
    message(paste0("k=", k, " was reduced to ", length(unique(test.labels)), " in test label prediction"))}
  
  return(list(test.labels = test.labels))
}

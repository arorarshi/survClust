
cv.voting<-function(fit,dd,kk, cmd.k=NULL, kk.test=TRUE, minlabel.test=TRUE){
  
  cc=combineDist(dd)
  if(is.null(cmd.k)){cmd.mat = cmdscale(cc, nrow(cc)-1, add=T)$points}
  if(!is.null(cmd.k)){cmd.mat = cmdscale(cc, cmd.k, add=T)$points}
  
  fit = fit[[kk-1]] 
  if(nrow(cmd.mat) != length(fit[[1]]$cv.labels)){stop("unequal samples in CV and cmd mat")}
  #remove fits with solutions != k
  ttt =unlist(lapply(fit, function(x) length(unique(x$cv.labels)) ))
  idx = which(ttt < kk)
  
  #test for solutions less than 5
  min.labels = lapply(fit, function(x) min(table(x$cv.labels)))
  idx2 = which(min.labels < 5)
  
  if(length(idx)==0){ kk.test=FALSE}
  if(length(idx2)==0){minlabel.test=FALSE}
  
  if(kk.test==TRUE & minlabel.test==TRUE){idx = unique(c(idx,idx2)); fit = fit[-idx]}
  if(kk.test==TRUE & minlabel.test==FALSE){fit = fit[-idx]}
  if(kk.test==FALSE & minlabel.test==TRUE){fit = fit[-idx2]}
  if(kk.test==FALSE & minlabel.test==FALSE){fit = fit}
  
  message(paste0("performing consensus on ", length(fit), " rounds"))
  if(length(fit)==0){stop("cross validation returned labels not equal to kk, pick another kk OR relax filters")}
  cv.rounds = length(fit)
  
  cmd.mat = cmd.mat[names(fit[[1]]$cv.labels),]
  
  centroids<-list()
  for (i in 1:cv.rounds){
    centroids[[i]]<-get.centroid(cmd.mat, fit[[i]]$cv.labels,i)
  }
  
  centroids.all<-do.call(rbind.data.frame, lapply(centroids, function(x) x))
  #do kmeans on the centroids
  centroids.kmeans<-kmeans(centroids.all,kk,nstart=100)
  all.cluster<-centroids.kmeans$cluster
  #print(table(all.cluster))
  #print(all.cluster)
  relabels<-list()
  for(i in 1:cv.rounds){
    pattern = paste0("f",i,"_k")
    relabels[[i]]<-get.relabel(pattern, fit[[i]]$cv.labels, all.cluster,kk)
    
  }
  relabels.all<-do.call(rbind.data.frame, lapply(relabels, function(x) x))
  relabels.all = apply(relabels.all, 2, as.numeric)
  colnames(relabels.all) = names(fit[[1]]$cv.labels)
  final.labels<-apply(relabels.all,2,function(x) names(table(x))[which.max(table(x))])
  return(unlist(final.labels))
}


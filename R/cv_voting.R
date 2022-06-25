
#' For a \code{survClust} fit, return consolidated labels across rounds of cross validation for a specific \code{k}. 
#' Note that cv.fit already has consolidated class labels across folds
#'
#' @param cv.fit fit objects as returned from \code{cv.survclust} 
#' @param dat.dist weighted distance matrices from \code{getDist}
#' @param pick_k choice of k cluster to summarize over rounds of cross validation 
#' @param cmd.k number of dimensions used by \code{cmdscale} to perform clustering on. Defaults is \code{n-1}  
#' @param pick_k.test logical, only selects cv.fit solutions where the resulting solution after consolidation contains \code{pick_k} classes. Default TRUE. Avoids edge cases, but in some cases FALSE might be desirable  
#' @param minlabel.test logical, only selects cv.fit solutions where classes have a minimum of 5 samples. Default TRUE. Avoids edge cases, but in some cases FALSE might be desirable
#'
#' @return final.labels consolidated class labels over rounds of cross-validation 
#' @export
#'
#' @examples

cv_voting<-function(cv.fit,dat.dist,pick_k, cmd.k=NULL, pick_k.test=TRUE, minlabel.test=TRUE){
  
  cc=combineDist(dat.dist)
  if(is.null(cmd.k)){cmd.mat = cmdscale(cc, nrow(cc)-1, add=T)$points}
  if(!is.null(cmd.k)){cmd.mat = cmdscale(cc, cmd.k, add=T)$points}
  
  cv.fit = cv.fit[[pick_k-1]] 
  if(nrow(cmd.mat) != length(cv.fit[[1]]$cv.labels)){stop("unequal samples in CV and cmd mat")}
  #remove cv.fits with solutions != k
  ttt =unlist(lapply(cv.fit, function(x) length(unique(x$cv.labels)) ))
  idx = which(ttt < pick_k)
  
  #test for solutions less than 5
  min.labels = lapply(cv.fit, function(x) min(table(x$cv.labels)))
  idx2 = which(min.labels < 5)
  
  if(length(idx)==0){ pick_k.test=FALSE}
  if(length(idx2)==0){minlabel.test=FALSE}
  
  if(pick_k.test==TRUE & minlabel.test==TRUE){idx = unique(c(idx,idx2)); cv.fit = cv.fit[-idx]}
  if(pick_k.test==TRUE & minlabel.test==FALSE){cv.fit = cv.fit[-idx]}
  if(pick_k.test==FALSE & minlabel.test==TRUE){cv.fit = cv.fit[-idx2]}
  if(pick_k.test==FALSE & minlabel.test==FALSE){cv.fit = cv.fit}
  
  message(paste0("performing consensus on ", length(cv.fit), " rounds"))
  if(length(cv.fit)==0){stop("cross validation returned labels not equal to pick_k, pick another pick_k OR relax filters")}
  cv.rounds = length(cv.fit)
  
  cmd.mat = cmd.mat[names(cv.fit[[1]]$cv.labels),]
  
  centroids<-list()
  for (i in 1:cv.rounds){
    centroids[[i]] <- .get_centroid(cmd.mat, cv.fit[[i]]$cv.labels,i)
  }
  
  centroids.all<-do.call(rbind.data.frame, lapply(centroids, function(x) x))
  #do kmeans on the centroids
  centroids.kmeans<-kmeans(centroids.all,pick_k,nstart=100)
  all.cluster<-centroids.kmeans$cluster
  #print(table(all.cluster))
  #print(all.cluster)
  relabels<-list()
  for(i in 1:cv.rounds){
    pattern = paste0("f",i,"_k")
    relabels[[i]] <- .get_relabel(pattern, cv.fit[[i]]$cv.labels, all.cluster,pick_k)
    
  }
  relabels.all<-do.call(rbind.data.frame, lapply(relabels, function(x) x))
  relabels.all = apply(relabels.all, 2, as.numeric)
  colnames(relabels.all) = names(cv.fit[[1]]$cv.labels)
  final.labels<-apply(relabels.all,2,function(x) names(table(x))[which.max(table(x))])
  return(unlist(final.labels))
}


#' Calculate standardized pooled within-cluster sum of squares (SPWSS)
#'
#' @param dist_mat a matrix of distance/dissimilarities on which to compute SPWSS
#' @param labels class labels 
#'
#' @return spwss metric for specified class labels over distance matrix
#' @export

get_spwss_stats<-function(dist_mat,labels){
  ll = unique(labels)
  dist_mat[lower.tri(dist_mat,diag=TRUE)]<-NA
  tss = sum(dist_mat, na.rm=T)
  wss<-rep(NA, length(ll))
  for (i in 1:length(ll)){
    wss[i] = sum(dist_mat[labels==ll[i], labels==ll[i]], na.rm=T)
  }
  
  tot.withinss = sum(wss)
  within.over.tot = tot.withinss/tss
  return((spwss = within.over.tot))
}

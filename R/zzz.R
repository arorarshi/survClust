
#' survClust
#' 
#' survClust is an outcome weighted integrative clustering algorithm used to classify multi-omic samples on their available time to event information. 
#' 
#' @docType package
#' @import Rcpp MultiAssayExperiment pdist survival
#' @importFrom Rcpp evalCpp
#' @importFrom pdist pdist
#' @importFrom survival survdiff
#' @importFrom survival survfit 
#' @useDynLib survClust
#' @name survClust
NULL 
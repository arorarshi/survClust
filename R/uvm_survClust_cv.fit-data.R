#' survClust cv.survclust output of integrated TCGA UVM Mutation and Copy Number datasets. 
#' 
#' 
#' The output is a list object consisting of 6 sub-lists for \code{k = 2:7}, with 10 \code{cv.survclust} outputs (for each round of cross-validation), each consisting of \code{cv.labels, cv.logrank, cv.spwss} for 3 folds. 
#' 
#' @docType data
#'
#' @usage data(uvm_survClust_cv.fit)
#'
#' @format An object of class \code{"list"}
#'
#' @keywords datasets
#'
#' @examples
#' data(uvm_survClust_cv.fit)
#' names(uvm_survClust_cv.fit[[1]][[1]])
"uvm_survClust_cv.fit"
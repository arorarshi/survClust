#' TCGA UVM Mutation and Copy Number datasets 
#'
#' A list of length 2 with 
#' TCGA UVM Mutation data with 80 samples and 87 genes
#' TCGA UVM Copy Number data with 80 samples and 749 segments. See Appendix
#' in vignette for more details. Teh data is downloaded from here
#' https://gdc.cancer.gov/about-data/publications/pancanatlas 
#' 
#' @docType data
#'
#' @usage data(uvm_dat)
#'
#' @format An object of class \code{"list"}
#'
#' @keywords datasets
#'
#' @examples
#' data(uvm_dat)
#' uvm_dat[[1]][1:5,1:5]
#' uvm_dat[[2]][1:5,1:5]
"uvm_dat"
#' Simulated survival dataset with accompanying \code{simdat}
#'
#' A matrix with simulate time-event data with 150 samples x 2 columns with a 
#' 3-class structure with median survival of 4.5, 3.25 and 2 yrs respectively.  
#' such that 15 features are distinct and associated with survival, other 15 
#' features are just distinct and not associated with survival and remaining 120
#' are noise. See how this dataset was generated in the vignette
#'
#' @docType data
#'
#' @usage data(simsurvdat)
#'
#' @format An object of class \code{"matrix"}
#'
#' @keywords datasets
#'
#' @examples
#' data(simsurvdat)
#' dim(simsurvdat)
#' head(simsurvdat)
"simsurvdat"
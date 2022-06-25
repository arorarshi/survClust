
#' returns cross-validated labels for optimal \code{pick_k}
#'
#' @description 
#' After plotting boxplots via \code{plotStats}, choose the optimal \code{pick_k} and pass it through \code{consensus_summary} 
#' to get consensus class labels from cross validation via \code{cv.survclust} 
#' 
#' @param pick_k choice of optimal \code{k} from \code{plotStats}
#' @param dist.dat pass the output of \code{getDist}
#' @param cv.fit output of \code{cv.survclust}
#' @param survdat A matrix, containing two columns - 1st column \code{time} and 2nd column containing \code{events} information.
#' @param plot.km logical, to plot Kaplan Meier curve of the cross validated labels of optimal k
#' @param ... Additional parameters as passed to \code{boxplot}
#'
#' @return
#' \itemize{
#'  \item{cv.labels}{cross validated labels after consolidating cross validation results of optimal \code{k}}
#'  }
#'  
#' @author Arshi Arora
#' 
#' @export
#' 
consensus_summary<-function(pick_k, dist.dat, cv.fit, survdat, plot.km=TRUE,...){
  
  cv.labels = cv_voting(cv.fit, dist.dat, pick_k)
  print("cross validated labels:")
  print(table(cv.labels))
  
  if(plot.km==TRUE){
    
    inter<-intersect(rownames(survdat), names(cv.labels))
    ff = cv.labels[inter]; survdat = survdat[inter,]
    ff = as.factor(ff)
    km=survfit(Surv(survdat[,1], survdat[,2])~ff)
    plot(km,mark.time=T,xlab="Time",col=1:length(unique(cv.labels)),main="cross validated solution",bty="l", ...)
    legend("bottomleft", paste0("c",1:length(unique(cv.labels))), bty="n", lty=1, col=1:length(unique(cv.labels)))
    
  }
  
  return(cv.labels)
}
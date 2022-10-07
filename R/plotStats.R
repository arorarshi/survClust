

#' Plot the output from \code{getStats}
#'
#' @param out.getStats list output from \code{getStats} 
#' @param labels labels to print on the boxplot. Default is \code{2:8} 
#' @param ... {additional arguments as passed to \code{boxplot} function}
#'
#' @details plots boxplots summarizing output of \code{cv.survclust} calculated via \code{getStats}. 
#' Use this to pick optimal \code{k}. Optimal \code{k} maximized logrank and minimizes SPWSS similar to the elbow 
#' method. Use \code{consensus_summary} to pick the best \code{k} and arrive at unique consolidated class labels
#' @return a plot with three boxplots summarizing logrank, standardized pooled within sum of squares (SPWSS) 
#' and if any class label has less than 5 samples
#' 
#' @examples 
#' library(survClust)
#' ss_stats <- getStats(uvm_survClust_cv.fit, kk=7, cvr=10)
#' plotStats(ss_stats, 2:7)
#' 
#' @author Arshi Arora
#' @export
plotStats <- function(out.getStats,labels=NULL,...){
    
    if(is.null(labels)){labels <- 2:8}
    dots <- list(...)
    #if(is.null(dots$main)){dots$main = paste0(nrow(out.getStats$lr)," datapoints") }
    
    par(mfrow=c(2,2))
    boxplot(out.getStats$lr, frame.plot=FALSE, names=labels, ylab="logrank",outline=FALSE,...)
    
    #boxplot of within over tot across datasets
    boxplot(out.getStats$spwss, frame.plot=FALSE, names=labels, ylab="SPWSS", ...)
    
    plot(out.getStats$bad.sol, frame.plot=FALSE,type='o', pch=8, ylab="solutions <= 5",xlab="k",xaxt="n",lwd=2,cex=1,...)
    axis(1, at=1:length(labels), labels=labels)
}




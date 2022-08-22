
#' Compute fit statistics after cross validation via \code{cv.survclust}
#'
#' @param cv.fit output from \code{cv.survclust} object 
#' @param kk number of \code{k} clusters on which \code{cv.survclust} was run, default is \code{8}
#' @param cvr round of cross-validation on which \code{cv.survclust}  was run, default is \code{50} 
#' @details \code{getStats} calculates Logrank statistic and standardized pooled within sum of squares (SPWSS) across 
#' cross-validated labels. Visualize it via \code{plotStats}
#' @return A list of the following
#' \itemize{
#'  \item{lr} {log rank statistic}
#'  \item{spwss} {standardized pooled within sum of squares}
#'  \item{bad.sol} {number of solutions for each \code{kk} that have cluster class \code{<5} samples}
#'  }
#'  @author Arshi Arora
#' @export
getStats<-function(cv.fit,kk=8, cvr=50){
  kk <- kk-1
  cvr <- cvr
  lr <- matrix(NA, nrow=cvr, ncol=kk)
  spwss <- matrix(NA, nrow=cvr, ncol=kk)
  bad.sol <- rep(NA, length=kk)
  
  for(i in 1:kk){
    lr[,i] <- t(unlist(lapply(cv.fit[[i]], function(x) x$cv.logrank)))
    spwss[,i] <- t(unlist(lapply(cv.fit[[i]], function(x) x$cv.spwss)))
    min.labels <- unlist(lapply(cv.fit[[i]], function(x) min(table(x$cv.labels)) [1]))
    idx <- which(min.labels <=5)
    #if (length(idx) !=0){llk = llk[idx]
    #    print(paste0("k= ", i, " has few samples in cluster solution for CV round - ",idx))}
    bad.sol[i]<-length(idx)
    
  }
  return(list(lr=lr, spwss = spwss, bad.sol = bad.sol))
}


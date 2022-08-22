
#' Integrates weighted distance matrices
#'
#'@description
#'\code{combineDist} integrates weighted distances matrices from \code{getDist}. 
#'All data types are now collapsed into one \code{NxN} matrix. 

#' @param dist.dat  list of weighted data matrices from \code{getDist}
#' @details 
#' \code{combineDist} integrates and does cleaning of missing pair of samples. 
#' if \code{datasets} list had non-overlapping samples, 
#' then \code{combineDist} retains only those samples that have full information after accounting for all data types.
#'
#' @return
#' \itemize{
#' \item{combMatFull}{A matrix. Combine normalized information across \code{m} genomic data types into \code{NxN} matrix, 
#' where \code{N} is the union of all samples across \code{m} data types/ or samples with complete pairwise information. 
#' Final matrix should not have any NAs}
#' }
#' @author Arshi Arora
#' 
#'
#' @export
#' 
combineDist<-function(dist.dat){
    
    arr <- array(unlist(dist.dat), dim = c(nrow(dist.dat[[1]]),ncol(dist.dat[[1]]),length(dist.dat)))
    #calculate mean distance after removing NA
    combMat <- rowMeans(arr, dim = 2, na.rm = TRUE)
    rownames(combMat) <- rownames(dist.dat[[1]])
    colnames(combMat) <- colnames(dist.dat[[1]])
    #copy the rest of the triangle
    combMat[lower.tri(combMat)] <- t(combMat)[lower.tri(combMat)]
    colna <- colSums(is.na(combMat))
    tab.idx <- sort(table(unlist(apply(combMat,2,function(x) which(is.na(x))))),decreasing=TRUE)
    idx <- as.numeric(names(tab.idx))
    iter <- 1
    combMatFull <- combMat
    #remove incomplete pairwise information
    while(sum(is.na(combMatFull))!=0){
        if (iter==1){del.idx <- idx[iter]}
        if (iter!=1){del.idx <- c(del.idx, idx[iter])}
        combMatFull = combMat[-del.idx, -del.idx]
        iter <- iter + 1
    }
    
    if(nrow(combMatFull) < 2){stop("only one sample left after overlaping for complete pairs, low sample overlap!")}
    return(combMatFull)
}

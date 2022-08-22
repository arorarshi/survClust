

#' perform supervised clustering for a particular \code{k}
#'
#'@description
#'\code{survClust} function performs supervised clustering on a \code{combineDist} output for a particular \code{k}. 
#'It uses all \code{n-1} dimensions for clustering.

#' @param combine.dist integrated weighted distance matrix from \code{combineDist}
#' @param survdat A nx2 matrix consisting of survival data with \code{n} samples and first column as time and second column as events, with samples as rownames
#' @param k choice of \code{k} to perform clustering on samples
#' @param cmd.k number of dimensions used by \code{cmdscale} to perform clustering on samples. Defaults is \code{n-1} 
#'
#' @return
#' \itemize{
#'  \item{fit} { returns a list , \code{fit} consisting of all clustering samples as in \code{kmeans} 
#'  \code{fit.lr}, computed logrank statistic between \code{k} clusters}
#' }
#' 
#' @author Arshi Arora
#' @export
survClust<-function(combine.dist,survdat,k, cmd.k=NULL){
    if(is.null(rownames(survdat)))
        stop("rowanmes of survdat can't be NULL")
    
    if(is.null(rownames(combine.dist)))
        stop("rowanmes of combine.dist can't be NULL")
    
    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("survival package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    
    if (!requireNamespace("pdist", quietly = TRUE)) {
        stop("pdist package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    
    inter <- intersect(rownames(survdat), rownames(combine.dist))
    
    clin <- survdat[inter,]
    clin <- apply(clin,2,as.numeric)
    
    #run cmdscale
    combine.dist <- combine.dist[inter,inter]
    
    if(is.null(cmd.k)){cmd.k <- nrow(combine.dist)-1 }
    if(!(is.null(cmd.k))){cmd.k <- as.numeric(cmd.k) }
    
    cmd.combine.dist<-cmdscale(combine.dist,k=cmd.k, add=TRUE)$points
    
    my.k <- as.numeric(k)
    survobj <- Surv(clin[,1], clin[,2])
    
    #run kmeans with 100 starts
    fit <- kmeans(cmd.combine.dist,my.k,nstart=100)
    
    #caluclate logrank
    fit.lr <- survdiff(Surv(clin[,1], clin[,2]) ~ fit$cluster)$chisq
    fit$fit.lr <- fit.lr
    
    #return fit and its logrank
    return(fit)
}


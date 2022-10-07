
#Compute predicted test labels on survClust fit object 
.predict_test_label <- function(all.cmd,fit,k){
    all.cmd <- as.matrix(all.cmd)
    train.snames <- names(fit$cluster)
    test.snames <- setdiff(rownames(all.cmd),train.snames)
    
    #where row - samples, col - genes
    centroid <- matrix(NA, nrow = k, ncol = ncol(all.cmd))
    for (kk in 1:k) {
        #meaning k clust has one sample. #WARNING #check
        if(is.vector(all.cmd[names(fit$cluster)[which(fit$cluster==kk)],]) & ncol(all.cmd) > 1){
            message("k=",k, " training cluster has one sample, prediction might be inaccurate")
            centroid[kk, ] <- all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ]
        }
        
        if (!(is.null(dim(all.cmd[names(fit$cluster)[fit$cluster==kk],])))){
            if(ncol(all.cmd)> 1){centroid[kk, ] <- apply(all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ], 2, mean)}
        }
        
        if(ncol(all.cmd)==1){centroid[kk,] <- mean(all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ])}
    }
    
    dist.whole <- apply(centroid,1,function(x) as.matrix(pdist::pdist(x,all.cmd)))
    
    #assign the cluster membership
    dist.labels <- apply(dist.whole,1,which.min)
    names(dist.labels) <- rownames(all.cmd)
    test.labels <- dist.labels[test.snames]
    
    #is missing a class label via pdist
    if(length(unique(test.labels)) != k){
        message("k=", k, " was reduced to ", length(unique(test.labels)), " in test label prediction")}
    
    return(list(test.labels = test.labels))
}

.cv_relabel<-function(mat, train.labels,cv.test.labels, k,fold){
    
    centroids<-list()
    
    for(i in seq_along(train.labels)){
        centroids[[i]] <- .get_centroid(mat, train.labels[[i]],i)
    }
    
    centroids.all <- do.call(rbind.data.frame, lapply(centroids, function(x) x))
    #do kmeans on the centroids
    centroids.kmeans <- kmeans(centroids.all,k,nstart=20)
    #print(centroids.kmeans$cluster)
    #centroids cluster labels
    all.cluster <- centroids.kmeans$cluster
    
    relabel <- rep(NA,nrow(mat))
    names(relabel) <- rownames(mat)
    
    for(i in seq_len(fold)){
        pattern <- paste0("f",i,"_k")
        rr <- .get_relabel(pattern, cv.test.labels[[i]], all.cluster,k)
        relabel[names(rr)] <- rr
    }
    
    return(relabel)
}

#' performs cross validation on supervised clustering, \code{survClust} for a particular \code{k}. \code{cv.survclust} runs 
#' 
#' @description 
#'\code{cv.survclust} performs \code{k} fold cross-validation, runs \code{survClust} on each training and 
#'hold out test fold and return cross-validated supervised cluster labels.  
#'
#' @param datasets A list object containing \code{m} data matrices representing \code{m} different genomic data types measured in a set of \code{N~m} samples. 
#' OR \code{\link{MultiAssayExperiment}} object of desired types of data. 
#' For list of matrices, each matrix, the rows represent samples, and the columns represent genomic features. Each data matrix is allowed to have different samples
#' @param survdat A matrix, containing two columns - 1st column \code{time} and 2nd column containing \code{events} information.
#' OR this information can be provided as a part of \code{colData} \code{MultiAssayExperiment}
#' @param k integer, choice of \code{k} to perform clustering on samples
#' @param fold integer,  number of folds to run cross validation 
#' @param cmd.k integer, number of dimensions used by \code{cmdscale} to perform clustering on samples. Defaults is \code{n-1}
#' @param type Specify \code{type="mut"}, if datasets is of length \code{1} and contains \code{binary} data only.
#'
#' @return
#' \itemize{
#'  \item{cv.labels}{returns cross validated class labels for \code{k} cluster} 
#'  \item{cv.logrank}{logrank test statistic of cross validated label}  
#'  \item{cv.spwss}{standardized pooled within-cluster sum of squares calculated from cross-validation  class labels }
#'  }
#'  
#' @examples 
#' library(survClust)
#' cv.fit <- cv.survclust(datasets = simdat, survdat = simsurvdat, k = 3, fold=3 )
#' 
#' @author Arshi Arora
#'
#' @export
cv.survclust <- function(datasets, survdat=NULL,k,fold, cmd.k=NULL, type=NULL){
    
    my.k <- as.numeric(k)
    fold <- as.numeric(fold)
    
    #checks for mae
    if(is(datasets, "MultiAssayExeriment")){survdat <- colData(datasets)}
    #mae forces survdat to have rownames. 
    
    if(is.null(survdat))
        stop("if datasets are provided as list of matrices, you need to provide survival data\n
           as a matrix OR MAE object lacks colData and survival information")
    if(length(unique(survdat[,2]))==1)
        stop("no deaths or censor events found in the survival data")
    
    # add other checks
    if (is.list(datasets) == TRUE & !(is(datasets,"MultiAssayExperiment"))){
        #convert everything to numeric - force user to provide a numeric matrix
        #dat<-lapply(datasets, function(x) as.data.frame(aaply(x,1,as.numeric,.drop=FALSE)) )
        rnames <- unlist(lapply(datasets, function(x) rownames(x)))
        rnames <- unique(rnames)
        
        if(is.null(rnames))
            stop("rowanmes=NULL, add sample names to matrix of datasets list object")
        
        if(is.null(rownames(survdat)))
            stop("rowanmes(survdat) cannot be NULL")
    }
    
    dist.dat <- getDist(datasets,survdat, type=type)
    #this for calculating ss on test labels
    combine.dist <- combineDist(dist.dat)
    inter <- intersect(rownames(survdat), rownames(combine.dist))
    ll <- seq(1,length(inter))
    
    #we will use this for cluster relabeling of test
    if(is.null(cmd.k)){this.k <- nrow(combine.dist)-1}
    if(!(is.null(cmd.k))){this.k <- as.numeric(cmd.k)}
    
    combine.dist.cmd <- cmdscale(combine.dist, k=this.k, add=TRUE)$points
    clin <- survdat[inter,]
    clin <- apply(clin,2,as.numeric)
    clin.whole <- clin
    rownames(clin.whole) <- inter
    
    folds <- cut(seq(1,length(ll)),breaks=fold,labels=FALSE)
    ll.rand <- sample(ll,length(ll))
    
    #Perform n fold cross validation
    cv.test.labels <- list()
    #cv.test.rand.index =NA
    survclust_fit <- list()
    
    for(i in seq_len(fold)){
        #Segement your data by fold using the which() function
        test.idx <- ll.rand[which(folds==i)]
        train.idx <- setdiff(ll,test.idx)
        train.snames <- rownames(clin.whole)[train.idx]
        clin.train <- clin.whole[train.snames,]
        
        #multiply by coxph abs(log(HR))
        distwt <- getDist(datasets,survdat,cv=TRUE,train.snames, type=type)
        train.dist.dat <- distwt$train
        all.dist.dat <- distwt$all
        
        #combine dist
        train.combine.dist <- combineDist(train.dist.dat)
        all.combine.dist <- combineDist(all.dist.dat)
        inter <- intersect(rownames(survdat), rownames(all.combine.dist))
        all.combine.dist <- all.combine.dist[inter,inter]
        
        if(is.null(cmd.k)){cmd.k.all <- nrow(all.combine.dist)-1 }
        if(!(is.null(cmd.k))){cmd.k.all <- as.numeric(cmd.k) }
        #as cmd is on dist, and nrow is different for all and training set
        #but multiplied with training HR
        cmd.whole <- cmdscale(all.combine.dist,k=cmd.k.all, add=TRUE)$points
        #get training fit labels
        fit <- survClust(train.combine.dist,survdat,my.k,cmd.k)
        
        #calculate test logrank and concordance
        #we basically predict on whole
        test <- .predict_test_label(cmd.whole,fit,my.k)
        cv.test.labels[[i]] <-  test$test.labels
        survclust_fit[[i]] <- fit
    }
    
    train.labels <- lapply(survclust_fit, function(x) x$cluster)
    cv.test.relabels <- .cv_relabel(combine.dist.cmd, train.labels,cv.test.labels,my.k,fold)
    min.labels <- min(table(cv.test.relabels))
    idx <- which(min.labels <=5)
    if (length(idx)!=0){message("k= ", my.k, " has 5 or few samples in cluster solution")}
    
    message("finished ", fold, " cross validation, total samples-", length(cv.test.relabels))
    cv.test.relabels <- cv.test.relabels[rownames(clin.whole)]
    if (length(unique(cv.test.relabels)) != my.k){warning("Test labels not equal to chosen k ",my.k) }
    
    #if everything collapses after test relabeling
    if (length(unique(cv.test.relabels)) ==1){
        cv.all.logrank <- NA
        warning("Everything collapsed after test relabel, logrank test is NA")
    }
    
    if (length(unique(cv.test.relabels)) !=1){cv.all.logrank <- survival::survdiff(Surv(clin.whole[names(cv.test.relabels),1], clin.whole[names(cv.test.relabels),2]) ~ cv.test.relabels)$chisq}
    #cv.all.conc <- summary(coxph(Surv(clin.whole[names(cv.test.relabels),1], clin.whole[names(cv.test.relabels),2]) ~ cv.test.relabels))$concordance[1]
    cv.test.spwss <- get_spwss_stats(combine.dist, cv.test.relabels)
    cv.fit <- list(cv.labels = cv.test.relabels, cv.logrank = cv.all.logrank, cv.spwss = cv.test.spwss)
    return(cv.fit)
    
}


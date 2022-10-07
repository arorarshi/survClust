
.getWeights<-function(mat, survdat,cv=FALSE, train.snames=NULL){
    #how can we compute faster HRs, apply is faster than for 
    
    if(!cv){
        
        #always intersect with survdat
        inter <- intersect(rownames(survdat),rownames(mat))
        mat <- mat[inter,,drop=FALSE]
        survdat <- survdat[inter,]
        
        #calculate survobj upfront
        #set infinite HR to a small value, and NA to 0
        survobj <- survival::Surv(as.numeric(survdat[,1]),as.numeric(survdat[,2]))
        logHR <- apply(mat,2, function(y) tryCatch({
            summary(survival::coxph(survobj ~ y))$coefficients[1]
        }, warning = function(w) {
            if(grepl("infinite", w$message))
            {message(w$message, " setting weights=0");return(1e-6)}
            
        }, error = function(e) {
            if(grepl("NA/NaN/Inf", e$message))
            {message(e$message, " .Setting weights=0"); return(0)}
            
        }) )
        logHR <- abs(logHR)
        #substitute it as 0. mostly should not happen, if it does meaning there is a gene that is not mutated  
        logHR[which(is.na(logHR))] <- 0
        
        mat.wt <- t(t(mat) * sqrt(logHR))
        return(mat.wt)
    }
    
    if(cv){
        inter <- intersect(rownames(survdat), intersect(rownames(mat),train.snames))
        mat.train <- mat[inter,,drop=FALSE]
        survdat <- survdat[inter,]
        
        #calculate survobj upfront
        survobj <- survival::Surv(as.numeric(survdat[,1]),as.numeric(survdat[,2]))
        logHR <- apply(mat.train,2, function(y) tryCatch({
            summary(survival::coxph(survobj ~ y))$coefficients[1]
        }, warning = function(w) {
            if(grepl("infinite", w$message))
            {message(w$message, " setting weights=0");return(0)}
            
        }, error = function(e) {
            if(grepl("NA/NaN/Inf", e$message))
            {message(e$message, " .Setting weights=0"); return(0)}
            
        }) )
        logHR <- abs(logHR)
        #NA in HR can result when a gene is not mutated
        #substitute it as 0 
        logHR[which(is.na(logHR))] <- 0 
        #on features which are columns, same in all and train
        mat.train.wt <- t(t(mat.train) * sqrt(logHR))
        #multiply HR on whole
        all.wt <- t(t(mat) * sqrt(logHR))
        mat.wt <- list(all=all.wt, train=mat.train.wt)
        return(mat.wt)
    }
    
}

.getWeights_mae <- function(datasets, cv, train.snames){
    
    mae_assays <- assays(datasets)
    dat_wt_mae <- lapply(mae_assays, 
                         function(x) .getWeights(mat=t(x), survdat = colData(datasets), cv, train.snames))
    return(dat_wt_mae)
    
}

.getUnionDist <- function(rnames,dat, type=NULL){
    
    dist.dat<-list()
    
    for (i in 1:length(dat)){
        m <- dat[[i]][intersect(rownames(dat[[i]]),rnames),,drop=FALSE]
        #if there are missing samples, add row of NAs for it
        if (length(intersect(rownames(m),rnames)) != length(rnames)){
            m.na <- matrix(NA,nrow=length(setdiff(rnames,rownames(m))),ncol=ncol(m))
            rownames(m.na) <- setdiff(rnames,rownames(m))
            m <- rbind(m,m.na)
            m <- m[rnames,]
        }
        if(!(is.null(type))){
            if(type=="mut"){
                dist.dat[[i]] <- dist_wtbinary(m)
                rownames(dist.dat[[i]]) = colnames(dist.dat[[i]]) = rownames(m)
            }
        }
        
        if(is.null(type)){    
            m2 <- m*m
            m2.ss <- sum(m2, na.rm=TRUE)
            m.tr <- m/sqrt(m2.ss)
            dist.dat[[i]] <- as.matrix(dist(m.tr, method="euclidean"))
        }
    }
    return(dist.dat)
}



.get_relabel <- function(pattern, olabel, centroid.cluster,kk){
    relabel <- rep(NA, length(olabel))
    names(relabel) <- names(olabel)
    
    for(i in 1:kk){
        kpattern <- paste0(pattern,i)
        idx <- which(names(centroid.cluster) == kpattern)
        #change current label to this
        if(length(idx)!=0){
            change.label <- centroid.cluster[idx]
            idx2 <- which(olabel == i)
            if(length(idx2)!=0){relabel[idx2] <- change.label}
        }
    }
    if(any(is.na(relabel))){warning("there is a NA in relabel, something is wrong with clustering, noisy data or pick lower k?")}
    return(relabel)
}


.get_centroid <- function (mat, labels, fold)
{
    ul <- unique(labels)
    if (ncol(mat) == 1) {
        warning("cmd reduces matrix to one eigen value! Noisy data?")
    }
    centroids <- matrix(NA, nrow = length(ul), ncol = ncol(mat))
    for (i in 1:length(ul)) {
        mat.ul <- mat[names(labels)[which(labels == ul[i])], ]
        if (is.vector(mat.ul)) {
            centroids[i, ] <- mat.ul
        }
        else {
            centroids[i, ] <- apply(mat.ul, 2, mean)
        }
    }
    rownames(centroids) <- paste0("f", fold, "_k", ul)
    return(centroids)
}


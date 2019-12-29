
###########################
# getDist
###########################

#computes weights for a given data, standalone
getWeights<-function(xx, survdat,cv=FALSE, train.snames=NULL){
  #add a check for missing rownames
  if(is.null(rownames(xx))){stop("provide rownames in getDist")}
  #how can we compute faster HRs, apply is faster than for 
  
  if(cv==FALSE){
    #always intersect with survdat
    inter = intersect(rownames(survdat),rownames(xx))
    xx = xx[inter,,drop=FALSE]
    survdat = survdat[inter,]
    
    #calculate survobj upfront
    survobj<-Surv(as.numeric(survdat[,1]),as.numeric(survdat[,2]))
    logHR<- apply(xx,2, function(y) tryCatch({
      summary(coxph(survobj ~ y))$coefficients[1]
    }, warning = function(w) {
      if(grepl("infinite", w$message))
      {message(paste0(w$message, " setting weights=0"));return(1e-6)}
      
    }, error = function(e) {
      if(grepl("NA/NaN/Inf", e$message))
      {message(paste0(e$message, " .Setting weights=0")); return(0)}
      
    }) )
    logHR = abs(logHR)
    #substitute it as 0. mostly should not happen, if it does meaning there is a gene that is not mutated  
    logHR[which(is.na(logHR))] = 0
    
    xx.wt <- t(t(xx) * sqrt(logHR))
    return(xx.wt)
  }
  
  if(cv==TRUE){
    inter = intersect(rownames(survdat), intersect(rownames(xx),train.snames))
    xx.train = xx[inter,,drop=FALSE]
    survdat = survdat[inter,]
    
    #calculate survobj upfront
    survobj<-Surv(as.numeric(survdat[,1]),as.numeric(survdat[,2]))
    logHR<- apply(xx.train,2, function(y) tryCatch({
      summary(coxph(survobj ~ y))$coefficients[1]
    }, warning = function(w) {
      if(grepl("infinite", w$message))
      {message(paste0(w$message, " setting weights=0"));return(0)}
      
    }, error = function(e) {
      if(grepl("NA/NaN/Inf", e$message))
      {message(paste0(e$message, " .Setting weights=0")); return(0)}
      
    }) )
    logHR = abs(logHR)
    #NA in HR can result when a gene is not mutated
    #substitute it as 0 
    logHR[which(is.na(logHR))] = 0 
    #on features which are columns, same in all and train
    xx.train.wt <- t(t(xx.train) * sqrt(logHR))
    #multiply HR on whole
    xx.wt <- t(t(xx) * sqrt(logHR))
    wt.mat<-list(all=xx.wt, train=xx.train.wt)
    return(wt.mat)
  }
  
}

getUnionDist<-function(rnames,dat, type=NULL){
  #compute distance across union of genes/features
  
  dist.dat<-list()
  
  for (i in 1:length(dat)){
    m <- dat[[i]][intersect(rownames(dat[[i]]),rnames),,drop=FALSE]
    #if there are missing samples, add row of NAs for it
    if (length(intersect(rownames(m),rnames)) != length(rnames)){
      m.na<-matrix(NA,nrow=length(setdiff(rnames,rownames(m))),ncol=ncol(m))
      rownames(m.na) = setdiff(rnames,rownames(m))
      m<-rbind(m,m.na)
      m<-m[rnames,]
    }
    if(!(is.null(type))){
      if(type=="mut"){dist.dat[[i]] = dist_wtbinary(m)
      rownames(dist.dat[[i]]) = colnames(dist.dat[[i]]) = rownames(m)}}
    
    if(is.null(type)){    
      m2 = m*m
      m2.ss = sum(m2, na.rm=T)
      m.tr = m/sqrt(m2.ss)
      dist.dat[[i]] = as.matrix(dist(m.tr, method="euclidean"))
    }
  }
  return(dist.dat)
}


getDist<-function(datasets,survdat,cv=FALSE,train.snames=NULL,type=NULL){
  # add other checks
  if (is.list(datasets) == F)
    stop("datasets must be a list")
  
  #if survdat has no events stop 
  if (length(unique(survdat[,2]))==1)
    stop("no deaths or censor events found")
  
  #convert everything to numeric
  dat<-lapply(datasets, function(x) as.data.frame(aaply(x,1,as.numeric,.drop=FALSE)) )
  rnames<-unlist(lapply(datasets, function(x) rownames(x)))
  rnames<-unique(rnames)
  
  if(is.null(rnames))
    stop("rowanmes=NULL, add sample names to matrix of datasets list object")
  
  dat.wt<-lapply(dat, function(x) getWeights(x,survdat,cv,train.snames))
  
  if(cv==TRUE){
    dat.all.wt<-lapply(dat.wt, function(x) x$all)
    dat.train.wt<-lapply(dat.wt, function(x) x$train)
    
    dat.all.dist<-getUnionDist(rnames, dat.all.wt,type)
    dat.train.dist<-getUnionDist(train.snames, dat.train.wt, type)
    dat.dist<-list(train=dat.train.dist, all=dat.all.dist)
    return(dat.dist)
  }
  
  if(cv==FALSE){
    dat.dist<-getUnionDist(rnames, dat.wt, type)
    return(dat.dist)
  }
  
}








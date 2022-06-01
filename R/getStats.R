
getStats<-function(cv.fit,kk=NULL, cvr=50){
  if(is.null(kk)){kk=7}
  cvr=cvr
  lr<-matrix(NA, nrow=cvr, ncol=kk)
  spwss<-matrix(NA, nrow=cvr, ncol=kk)
  bad.sol<-rep(NA, length=kk)
  
  for(i in 1:kk){
    lr[,i] = t(unlist(lapply(cv.fit[[i]], function(x) x$cv.logrank)))
    spwss[,i]<-t(unlist(lapply(cv.fit[[i]], function(x) x$cv.ss)))
    min.labels = unlist(lapply(cv.fit[[i]], function(x) min(table(x$cv.labels)) [1]))
    idx = which(min.labels <=5)
    #if (length(idx) !=0){llk = llk[idx]
    #    print(paste0("k= ", i, " has few samples in cluster solution for CV round - ",idx))}
    bad.sol[i]<-length(idx)
    
  }
  return(list(lr=lr, spwss = spwss, bad.sol = bad.sol))
}


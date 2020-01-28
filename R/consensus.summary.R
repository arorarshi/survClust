
getstats<-function(ll,kk=NULL, cvr=50){
  if(is.null(kk)){kk=7}
  cvr=cvr
  lr<-matrix(NA, nrow=cvr, ncol=kk)
  conc<-matrix(NA, nrow=cvr, ncol=kk)
  within.over.tot<-matrix(NA, nrow=cvr, ncol=kk)
  bad.sol<-rep(NA, length=kk)
  
  for(i in 1:kk){
    lr[,i] = t(unlist(lapply(ll[[i]], function(x) x$cv.logrank)))
    within.over.tot[,i]<-t(unlist(lapply(ll[[i]], function(x) x$cv.ss)))
    min.labels = unlist(lapply(ll[[i]], function(x) min(table(x$cv.labels)) [1]))
    idx = which(min.labels <=5)
    #if (length(idx) !=0){llk = llk[idx]
    #    print(paste0("k= ", i, " has few samples in cluster solution for CV round - ",idx))}
    bad.sol[i]<-length(idx)
    
  }
  return(list(lr=lr, within.over.tot=within.over.tot, bad.sol = bad.sol))
}


plotstats<-function(xx,labels=NULL,...){
  
  if(is.null(labels)){labels=2:8}
  dots<-list(...)
  if(is.null(dots$main)){dots$main = paste0(nrow(xx$lr)," datapoints") }
  
  par(mfrow=c(2,2))
  boxplot(xx$lr, frame.plot=F, names=labels, ylab="logrank",outline=F,...)
  
  #boxplot of within over tot across datasets
  boxplot(xx$within.over.tot, frame.plot=F, names=labels, ylab="SPWSS", ...)
  
  plot(xx$bad.sol, frame.plot=F,type='o', pch=8, ylab="solutions <= 5",xlab="k",xaxt="n",lwd=2,cex=1,...)
  axis(1, at=1:length(labels), labels=labels)
}


consensus.summary<-function(pick, dd, cv.fit, survdat, plot.km=TRUE,...){
  
  cv.labels = cv.voting(cv.fit, dd, pick)
  print("cross validated labels:")
  print(table(cv.labels))
  
  if(plot.km==TRUE){
    km.col = c("darkkhaki","darkmagenta","firebrick2","cornflowerblue","green","darkblue","black","hotpink2")
    
    inter<-intersect(rownames(survdat), names(cv.labels))
    ff = cv.labels[inter]; survdat = survdat[inter,]
    ff = as.factor(ff)
    km=survfit(Surv(survdat[,1], survdat[,2])~ff)
    plot(km,mark.time=T,xlab="Time",col=km.col[1:length(unique(cv.labels))],main="cross validated solution",bty="l", ...)
    legend("bottomleft", paste0("c",1:length(unique(cv.labels))), bty="n", lty=1, col=km.col[1:length(unique(cv.labels))])
    
  }
        
  return(cv.labels)
}


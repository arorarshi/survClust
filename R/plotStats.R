
plotStats<-function(ll,labels=NULL,...){
  
  if(is.null(labels)){labels=2:8}
  dots<-list(...)
  if(is.null(dots$main)){dots$main = paste0(nrow(xx$lr)," datapoints") }
  
  par(mfrow=c(2,2))
  boxplot(ll$lr, frame.plot=F, names=labels, ylab="logrank",outline=F,...)
  
  #boxplot of within over tot across datasets
  boxplot(ll$spwss, frame.plot=F, names=labels, ylab="SPWSS", ...)
  
  plot(ll$bad.sol, frame.plot=F,type='o', pch=8, ylab="solutions <= 5",xlab="k",xaxt="n",lwd=2,cex=1,...)
  axis(1, at=1:length(labels), labels=labels)
}




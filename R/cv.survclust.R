#Runs cross validation on survclust to attain best k
#calculate sum of squares of each simulated dataset
do.ss.stats<-function(mm,labels){
  ll = unique(labels)
  mm[lower.tri(mm,diag=TRUE)]<-NA
  tss = sum(mm, na.rm=T)
  wss<-rep(NA, length(ll))
  for (i in 1:length(ll)){
    wss[i] = sum(mm[labels==ll[i], labels==ll[i]], na.rm=T)
  }
  
  tot.withinss = sum(wss)
  within.over.tot = tot.withinss/tss
  return(within.over.tot)
}

get.centroid<-function (mat, labels, f)
{
  ul <- unique(labels)
  if (ncol(mat) == 1) {
    warning("cmd reduces matrix to one eigen value! Noisy data?")
  }
  centroids <- matrix(NA, nrow = length(ul), ncol = ncol(mat))
  for (i in 1:length(ul)) {
    mat.ul <- mat[names(labels)[which(labels == ul[i])], ]
    if (is.vector(mat.ul)) {
      centroids[i, ] = mat.ul
    }
    else {
      centroids[i, ] <- apply(mat.ul, 2, mean)
    }
  }
  rownames(centroids) = paste0("f", f, "_k", ul)
  return(centroids)
}

get.relabel<-function(pattern, olabel, centroid.cluster,kk){
  relabel<-rep(NA, length(olabel))
  names(relabel) = names(olabel)
  
  for(i in 1:kk){
    kpattern<-paste0(pattern,i)
    idx<-which(names(centroid.cluster)==kpattern)
    #change current label to this
    if(length(idx)!=0){
      change.label = centroid.cluster[idx]
      idx2 = which(olabel == i)
      if(length(idx2)!=0){relabel[idx2] = change.label}
    }
  }
  if(any(is.na(relabel))){warning("there is a NA in relabel, something is wrong with clustering, noisy data or pick lower k?")}
  return(relabel)
}

#this is done to provide meaning to cluster labels across folds.
#we run kmeans on centroid vector of all the folds to determine their closeness.
#random start -10
cv.relabel<-function(mat, train.labels,cv.test.labels, k,fold){
  
  centroids<-list()
  
  for(i in 1:length(train.labels)){
    centroids[[i]]<-get.centroid(mat, train.labels[[i]],i)
  }
  
  centroids.all<-do.call(rbind.data.frame, lapply(centroids, function(x) x))
  #do kmeans on the centroids
  centroids.kmeans<-kmeans(centroids.all,k,nstart=20)
  #print(centroids.kmeans$cluster)
  #centroids cluster labels
  all.cluster<-centroids.kmeans$cluster
  
  relabel<-rep(NA,nrow(mat))
  names(relabel) = rownames(mat)
  
  for(i in 1:fold){
    pattern = paste0("f",i,"_k")
    rr<-get.relabel(pattern, cv.test.labels[[i]], all.cluster,k)
    relabel[names(rr)] = rr
  }
  
  return(relabel)
}

#############################
# perform cross validation 
#############################

cv.survclust<-function(x, survdat,k,fold, cmd.k=NULL, type=NULL){
  
  my.k <- as.numeric(k)
  fold <- as.numeric(fold)
  
  #To get an idea of total samples
  dist.dat<-getDist(x,survdat, type=type)
  #this for calculating ss on test labels
  combine.dist<-combineDist(dist.dat)
  inter <- intersect(rownames(survdat), rownames(combine.dist))
  ll <- seq(1,length(inter))
  
  #we will use this for cluster relabeling of test
  if(is.null(cmd.k)){this.k = nrow(combine.dist)-1}
  if(!(is.null(cmd.k))){this.k = as.numeric(cmd.k)}
  
  combine.dist.cmd<-cmdscale(combine.dist, k=this.k, add=TRUE)$points
  clin<-survdat[inter,]
  clin<-apply(clin,2,as.numeric)
  clin.whole<-clin
  rownames(clin.whole) = inter
  
  folds <- cut(seq(1,length(ll)),breaks=fold,labels=FALSE)
  ll.rand<- sample(ll,length(ll))
  
  #Perform n fold cross validation
  cv.test.labels<-list()
  #cv.test.rand.index =NA
  survfit<-list()
  
  for(i in 1:fold){
    #Segement your data by fold using the which() function
    test.idx <- ll.rand[which(folds==i)]
    train.idx <- setdiff(ll,test.idx)
    train.snames = rownames(clin.whole)[train.idx]
    clin.train<- clin.whole[train.snames,]
    
    #multiply by coxph abs(log(HR))
    distwt<-getDist(x,survdat,cv=TRUE,train.snames, type=type)
    train.dist.dat<-distwt$train
    all.dist.dat<-distwt$all
    
    #combine dist
    train.combine.dist<-combineDist(train.dist.dat)
    all.combine.dist<-combineDist(all.dist.dat)
    inter <- intersect(rownames(survdat), rownames(all.combine.dist))
    all.combine.dist = all.combine.dist[inter,inter]
    
    if(is.null(cmd.k)){cmd.k.all =nrow(all.combine.dist)-1 }
    if(!(is.null(cmd.k))){cmd.k.all =as.numeric(cmd.k) }
    #as cmd is on dist, and nrow is different for all and training set
    #but multiplied with training HR
    cmd.whole = cmdscale(all.combine.dist,k=cmd.k.all, add=TRUE)$points
    #get training fit labels
    fit=survclust(train.combine.dist,survdat,my.k,cmd.k)
    
    #calculate test logrank and concordance
    #we basically predict on whole
    test =predict.test.label(cmd.whole,fit,my.k)
    cv.test.labels[[i]] = test$test.labels
    survfit[[i]]<-fit
  }
  
  train.labels = lapply(survfit, function(x) x$cluster)
  cv.test.relabels<-cv.relabel(combine.dist.cmd, train.labels,cv.test.labels,my.k,fold)
  min.labels = min(table(cv.test.relabels))
  idx = which(min.labels <=5)
  if (length(idx)!=0){message(paste0("k= ", my.k, " has 5 or few samples in cluster solution"))}
  
  message(paste0("finished ", fold, " cross validation, total samples-", length(cv.test.relabels)))
  cv.test.relabels = cv.test.relabels[rownames(clin.whole)]
  if (length(unique(cv.test.relabels)) != my.k){warning(paste0("Test labels not equal to chosen k ",my.k)) }
  
  #if everything collapses after test relabeling
  if (length(unique(cv.test.relabels)) ==1){
    cv.all.logrank = NA
    warning("Everything collapsed after test relabel, logrank test is NA")
  }
  
  if (length(unique(cv.test.relabels)) !=1){cv.all.logrank = survdiff(Surv(clin.whole[names(cv.test.relabels),1], clin.whole[names(cv.test.relabels),2]) ~ cv.test.relabels)$chisq}
  cv.all.conc = summary(coxph(Surv(clin.whole[names(cv.test.relabels),1], clin.whole[names(cv.test.relabels),2]) ~ cv.test.relabels))$concordance[1]
  cv.test.ss<-do.ss.stats(combine.dist, cv.test.relabels)
  cv.fit = list(cv.labels = cv.test.relabels, cv.logrank = cv.all.logrank, cv.concordance = cv.all.conc, cv.ss = cv.test.ss)
  return(cv.fit)
  
}


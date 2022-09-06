#' Calculates weighted distance matrix of multiple genomic data types 
#'
#' @description
#' Given multiple genomic data types (e.g., gene expression, copy number, DNA methylation, miRNA expression (continuous) and mutation (binary)) measured across samples, 
#' allowing for missing values (NA) and missing samples,  \code{getDist} calculates the survival weighted distance metric among samples. 
#' Used as an input to, \code{combineDist()}.
#' 
#' @param datasets A list object containing \code{m} data matrices representing \code{m} different genomic data types measured in a set of \code{N~m} samples. 
#' OR \code{\link{MultiAssayExperiment}} object of desired types of data. 
#' For list of matrices, each matrix, the rows represent samples, and the columns represent genomic features. Each data matrix is allowed to have different samples
#' @param survdat A matrix, containing two columns - 1st column \code{time} and 2nd column containing \code{events} information.
#' OR this information can be provided as a part of \code{colData} \code{MultiAssayExperiment}
#' @param cv logical. If \code{TRUE}, \code{train.names} cannot be NULL. Cross-validation will be performed on \code{train.names} samples, 
#' and the dataset will be split into training and test, and each respective matrices will be returned.
#' @param train.snames required if \code{cv=TRUE}. A vector of sample names treated as training samples.
#' @param type \code{NULL}. Specify \code{type="mut"}, if datasets is of length \code{1} and contains \code{binary} data only. See \code{details}
#'
#' @details
#' \code{getDist} allows for continuous and binary data type(s) in a matrix passed as a list. 
#' If the list only has a binary matrix data type. Set \code{type="mut"}. All data types are standardized internally. 
#' All data types are not expected to have common samples. Non-overlapping samples within data types are replaced with NA, and returned weighted matrix consists of union of all the samples.
#' @return
#' \itemize{
#'  \item{cv=FALSE,dist.dat}{returns a list of weighted data matrix/matrices, \code{dist.dat}}
#'  \item{cv=TRUE,dist.dat=list(train, all)}{ returns a list of training \code{train} weighted data matrix. 
#'  And the whole matrix weighed according to the weights computed on the training dataset \code{all}.  }
#'  }
#'  
#' @examples 
#' library(survClust)
#' dd <- getDist(simdat, simsurvdat)
#' 
#' @author Arshi Arora
#' 
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @export
getDist<-function(datasets,survdat=NULL,cv=FALSE,train.snames=NULL,type=NULL){
    
    if(is(datasets, "MultiAssayExperiment")){survdat <- colData(datasets)}
    if(is.null(survdat))
        stop("if datasets are provided as list of matrices, you need to provide survival data\n
           as a matrix OR MAE object lacks colData and survival information")
    if(length(unique(survdat[,2]))==1)
        stop("no deaths or censor events found in the survival data")
    
    # add other checks
    if (!(is(datasets,"MultiAssayExperiment"))){ if(!is.list(datasets)) stop("input data types in list")}
    if (is.list(datasets) & !(is(datasets,"MultiAssayExperiment"))){
        #convert everything to numeric - force user to provide a numeric matrix
        #dat<-lapply(datasets, function(x) as.data.frame(aaply(x,1,as.numeric,.drop=FALSE)) )
        rnames <- unlist(lapply(datasets, function(x) rownames(x)))
        rnames <- unique(rnames)
        
        if(is.null(rnames))
            stop("rowanmes=NULL, add sample names to matrix of datasets list object")
        
        if(is.null(rownames(survdat)))
            stop("rowanmes(survdat) cannot be NULL")
        
        dat.wt <- lapply(datasets, function(x) .getWeights(x,survdat,cv,train.snames))
    }
    
    
    
    #check if mae has been passed
    if(is(datasets,"MultiAssayExperiment")){
        
        dat.wt <- .getWeights_mae(datasets, cv, train.snames)
    }
    
    if(cv){
        dat.all.wt <- lapply(dat.wt, function(x) x$all)
        dat.train.wt <- lapply(dat.wt, function(x) x$train)
        
        dat.all.dist <- .getUnionDist(rnames, dat.all.wt,type)
        dat.train.dist <- .getUnionDist(train.snames, dat.train.wt, type)
        dat.dist <- list(train=dat.train.dist, all=dat.all.dist)
        return(dat.dist)
    }
    
    if(cv==FALSE){
        dat.dist <- .getUnionDist(rnames, dat.wt, type)
        return(dat.dist)
    }
    
    
    
}

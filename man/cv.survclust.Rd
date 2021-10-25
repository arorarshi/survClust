\name{cv.survclust}
\alias{cv.survclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
performs cross validation on supervised clustering, \code{survClust} for a particular \code{k}. \code{cv.survclust} runs 
}
\description{
\code{cv.survclust} performs \code{k} fold cross-validation, runs \code{survClust} on each training and hold out test fold and return cross-validated supervised cluster labels.  
}
\usage{
cv.survclust<-function(x, survdat,k,fold, cmd.k=NULL, type=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x} { A list object containing \code{m} data matrix/matrices representing \code{m} different genomic data types measured in a set of \code{N~m} samples. For each matrix, the rows represent samples, and the columns represent genomic features. Each data matrix is allowed to have different samples.}
\item{survdat} { A nx2 matrix consisting of survival data with \code{n} samples and first column as time and second column as events, with samples as rownames}
\item{k} { choice of \code{k} to cluster samples on}
\item{cmd.k} { number of dimensions used by \code{cmdscale} to perform clustering on. Defaults is \code{n-1}}
\item{type}{ \code{NULL}. Specify \code{type="mut"}, if datasets is of length \code{1} and contains \code{binary} data only.}
}

\value{
\item{cv.labels}{returns cross validated class labels for \code{k} cluster} \item{cv.logrank}{logrank test statistic of cross validated label}  
\item{cv.ss}{standardized pooled within-cluster sum of squares calculated from cross-validation  class labels }
}

\author{
Arshi Arora
}

\name{getStats}
\alias{getStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\code{getStats} computes necessary statistic after cross validation via \code{cv.survclust} 
}
\description{
\code{getStats} calculates log rank statistic and standardized pooled within sum of squares (SPWSS) across cross validated labels. Visualize it via \code{plotStats}
}
\usage{
getStats<-function(cv.fit,kk=NULL, cvr=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cv.fit}{output from \code{cv.survclust} object}
  \item{kk}{number of \code{k} clusters on which \code{cv.survclust} was run, default is \code{7}}
  \item{cvr}{ round of cross validation on which \code{cv.survclust} was run, default is \code{50}}
}
\value{
\item{lr}{log rank statistic}
\item{within.over.tot}{standardized pooled within sum of squares}
\item{bad.sol}{number of solutions for each \code{k} that have solutions \code{<5}}
}

\author{
Arshi Arora
}

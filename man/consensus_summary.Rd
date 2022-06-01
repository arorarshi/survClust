\name{consensus_summary}
\alias{consensus_summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
returns cross-validated labels for optimal \code{k}
}
\description{
After plotting boxplots via \code{plotStats}, choose the optimal \code{pick_k} and pass it through \code{consensus_summary} to get consensus class labels from cross validation via \code{cv.survclust} 
}
\usage{
consensus_summary(pick_k, dist.dat, cv.fit, survdat, plot.km=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pick_k}{choice of optimal \code{k} from \code{plotStats}}
  \item{dist.dat}{pass the output of \code{getDist}}
  \item{cv.fit}{output of \code{cv.survclust}}
  \item{survdat}{A matrix, containing two columns - 1st column \code{time} and 2nd column containing \code{events} information.}
  \item{plot.km}{logical, to plot Kaplan Meier curve of the cross validated labels of optimal k}
}

\value{
\item{cv.labels}{cross validated labels after consolidating cross validation results of optimal \code{k}}
}

\author{
Arshi Arora
}

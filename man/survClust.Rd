\name{survClust}
\alias{survClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{perform supervised clustering for a particular \code{k}
}

\description{
\code{survClust} function performs supervised clustering on a \code{combineDist} output for a particular \code{k}. It uses all \code{n-1} dimensions for clustering. 
}
\usage{
survClust(combine.dist,survdat,k, cmd.k=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{combine.dist} { integrated weighted distance matrix from \code{combineDist}}
\item{survdat} { A nx2 matrix consisting of survival data with \code{n} samples and first column as time and second column as events, with samples as rownames}
\item{k} { choice of \code{k} to cluster samples on}
\item{cmd.k} { number of dimensions used by \code{cmdscale} to perform clustering on. Defaults is \code{n-1}}
}

\value{
\item{fit} { returns a list , \code{fit} consisting of all clustering samples along with \code{fit.lr}, computed logrank statistic between \code{k} clusters}
}
\author{
Arshi Arora
}

\examples{

}

\name{plotstats}
\alias{plotstats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plots the output from \code{getstats}
}
\description{
plots boxplots summarizing output of \code{cv.survclust} calculated via \code{getstats}. Plots logrank, standardized within sum of squares and if any class label has less than 5 samples. 
}
\usage{
plotstats(ll, labels=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{ll}{output from \code{getstats}}
 \item{labels}{labels to print on the boxplot. Default is \code{2:8}}
 \item{...}{other arguments as passed to boxplot function}
}

\details{
Use \code{plotstats} to pick the optimal \code{k}. The one that maximizes logrank and minimized standardized pooled within sum of squares via the elbow method. Use \code{consensus.summary} to pick the optimal \code{k} and get one consolidated class labels. 
}

\author{
Arshi Arora
}

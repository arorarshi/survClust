\name{combineDist}
\alias{combineDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrates weighted distance matrices
}
\description{

\code{combineDist} integrates weighted distances matrices from \code{getDist}. All data types are now collapsed into one \code{NxN} matrix. 
}

\usage{
combineDist(dist.dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist.dat}{ list of weighted data matrices from \code{getDist}
}
}
\details{
\code{combineDist} integrates and does cleaning of missing pair of samples. if \code{datasets} list had non-overlapping smaples, then \code{combineDist} retains only those samples that have full information after accounting for all data types.

}
\value{
\item{combMatFull}{A matrix. Combine normalized information across \code{m} genomic data types into \code{NxN} matrix, where \code{N} is the union of all samples across \code{m} data types/ or samples with complete pairwise information. Final matrix should not have any NAs}

}
\author{
Arshi Arora
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

\dontrun{
dd = getDist(x,survdat)
cc = combineDist(dd)
}

}


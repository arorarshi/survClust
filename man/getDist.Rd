\name{getDist}
\alias{getDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates weighted distance matrix of multiple genomic data types 
}
\description{
\code{getDist} 

Given multiple genomic data types (e.g., gene expression, copy number, DNA methylation, miRNA expression and mutation(continuous and binary)) measured across samples, allowing for missing values (NA) and missing samples,  \code{getDist} calculates the survival weighted distance metric among samples. 

Used as an input to, \code{combineDist()}.
See also \code{survClust()}
}

\usage{
getDist(datasets, survdat, cv = FALSE, train.snames = NULL, type = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datasets}{ A list object containing \code{m} data matrices representing \code{m} different genomic data types measured in a set of \code{N~m} samples. For each matrix, the rows represent samples, and the columns represent genomic features. Each data matrix is allowed to have different samples. 
}
  \item{survdat}{ A matrix, containing two columns - 1st column \code{time} and 2nd column containing \code{events} information. 

}
  \item{cv}{ logical. If \code{TRUE}, \code{train.names} cannot be NULL. Cross-validation will be performed on \code{train.names} samples, and the dataset will be split into training and test, and each respective matrices will be returned. 

}
  \item{train.snames}{ required if \code{cv=TRUE}. A vector of sample names treated as training samples. 
}
  \item{type}{ \code{NULL}. Specify \code{type="mut"}, if datasets is of length \code{1} and contains \code{binary} data only. See \code{details}

}
}
\details{

\code{getDist} allows for continuous and binary data type(s) in a matrix passed as a list. If the list only has a binary matrix data type. Set \code{type="mut"}. All data types are standardized internally. 

All data types are not expected to have common samples. Non-overlapping samples within data types are replaced with NA, and returned weighted matrix consists of union of all the samples.
}

\value{

\item{cv=FALSE,dist.dat}{returns a list of weighted data matrix/matrices, \code{dist.dat}}

\item{cv=TRUE,dist.dat=list(train, all)}{ returns a list of training \code{train} weighted data matrix. And the whole matrix weighed according to the weights computed on the training dataset \code{all}.  }

}
\references{

1. Xing, E.P., Jordan, M.I., Russell, S.J. and Ng, A.Y., 2003. Distance metric learning with application to clustering with side-information. In Advances in neural information processing systems (pp. 521-528).

}
\author{
Arshi Arora
}
\note{
 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

######################
#simulate data type
######################

#simulation to test survival related + survival unrelated
#total samples = 150, total features = 150

set.seed(112)
n1 = 50 #class1
n2 = 50 #class2
n3 = 50 #class3
n = n1+n2+n3
p = 15 #survival related features (10%)
q = 120 #noise

#class1 ~ N(1.5,1), class2 ~ N(0,1), class3 ~ N(-1.5,1)

rnames = paste0("S",1:n)
x = NULL
x1a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
x2a = matrix(rnorm(n1*p), ncol=p)
x3a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[1]] = cbind(xa,xb)

## sample 15 other informant features, but scramble them.

permute.idx<-sample(1:length(rnames),length(rnames))
x1a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
x2a = matrix(rnorm(n1*p), ncol=p)
x3a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
xa = rbind(x1a,x2a,x3a)

x[[1]] = cbind(x[[1]],xa[permute.idx,])
rownames(x[[1]]) =  rnames

######################
#simulate survival data 
######################

set.seed(112)
l1 = log(2)/4.5
l2 = log(2)/3.25
l3 = log(2)/2
n = 50

s1<-rexp(n,rate=l1)
c1<-runif(n,0,10)
s2<-rexp(n,rate=l2)
c2<-runif(n,0,10)
s3<-rexp(n,rate=l3)
c3<-runif(n,0,10)

t1 = pmin(s1,c1)
t2 = pmin(s2,c2)
t3 = pmin(s3,c3)

event = c((t1==s1),(t2==s2),(t3==s3))
time = c(t1,t2,t3)
survdat = cbind(time, event)
rownames(survdat) = rnames
#plot(survfit(Surv(time,event) ~ truth), lty=1:3, mark.time=T, bty="l",lwd=1.5, main="simulated survival dataset")

dd = getDist(x,survdat)

}

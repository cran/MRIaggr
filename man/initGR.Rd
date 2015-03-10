\name{initGR}
\title{Growing Region initialization}
\alias{initGR}
\usage{
initGR(contrast,W,seed,range,range.seed,breaks,scale,trace,method)
}
\description{
Initialize arguments of the Growing Region algorithm. For internal use.
}
\arguments{
  \item{contrast}{the contrast value of each observation. \emph{numeric vector}.}
  \item{W}{the neighborhood matrix. \code{dgCMatrix}.}
  \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \code{positive integer vector} or \code{logical vector}.}
  \item{range}{the range of acceptable contrast values for the growing region group. \code{numeric vector of size 2}.}
  \item{range.seed}{the range of acceptable contrast values for the seeds. \code{numeric vector of size 2}.}
 
  \item{breaks}{the break points or the number of break points to use to categorize the contrast distribution. \code{numeric vector} or \code{postive integer}.}
  \item{scale}{should the contrast be scaled ? \code{logical}.}  
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}  
  
  \item{method}{the name of the function that called the initializer. \emph{character}.}
}
\concept{init.}
\keyword{functions,internal}
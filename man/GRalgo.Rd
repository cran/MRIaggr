\name{GRalgo}
\title{Growing Region algorithm}
\alias{GRalgo}
\description{
Perform the Growing Region algorithm proposed by (Revol et al., 1997). For internal use.
}
\usage{
GRalgo(contrast,W,seed,sigma_max,range,breaks,step,operator,iter_max,
         history_sigma,history_step,history_front)
}
\arguments{
  \item{contrast}{the contrast value of each observation. \emph{numeric vector}.}
  \item{W}{the neighborhood matrix. \code{dgCMatrix}.}
  \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \code{positive integer vector} or \code{logical vector}.}
  \item{sigma_max}{the maximum admissible value for the variability of the group contrast. \code{positive numeric}.}
  \item{range}{the range of acceptable contrast values for the growing region group. \code{numeric vector of size 2}.}
 
  \item{breaks}{the break points to use to categorize the contrast distribution. \code{numeric vector}.}
  \item{step}{the step between two consecutive breaks. \emph{numeric}.}  
  \item{iter_max}{the maximum number of iterations for the expansion of the growing region. \code{postive integer}.}  
  \item{operator}{should the median absolute deviation be used to estimte the variability of the group contrast (\code{"mad"}) or the standard deviation (\code{"sd"}).}  
  
  \item{history_sigma}{should the values of sigma be recorded ? \emph{logical}.}
  \item{history_step}{should the number of observations included in the GR set be recorded ? \emph{logical}.}
  \item{history_front}{should the propagation front of the GR set be recorded ? \emph{logical}.}
}
\references{
Chantal Revol and Michel Jourlin. \emph{A new minimum varance region growing algorithm for image segmentation}. Pattern Recognition Letters, 18(3):249-258,1997.
}
\concept{calc.}
\keyword{functions,internal}
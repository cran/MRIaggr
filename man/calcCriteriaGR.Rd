\name{calcCriteriaGR}
\title{Assessment of clustering quality}
\alias{calcCriteriaGR}
\description{
Compute several quality indexes of a two group clustering. For internal use.
}
\usage{
calcCriteriaGR(contrast,groups,W=NULL,sigma=NULL,breaks,criterion_d1=FALSE,
         criterion_entropy=TRUE,criterion_Kalinsky=TRUE,criterion_Laboure=TRUE)
}
\arguments{
  \item{contrast}{the contrast value of each observations. \emph{numeric vector}. REQUIRED.}
  \item{groups}{the indicator of group membership. \code{logical vector}. REQUIRED.}
  \item{W}{the neighborhood matrix. \code{dgCMatrix} or \code{NULL} leading to not compute the d1 criterion.}
  \item{sigma}{the \code{sigma_max} that have been used in the GR algorithm. \code{positive numeric vector}.}
  \item{breaks}{the break points to use to categorize the contrast distribution. \code{numeric vector}.}
 
  \item{criterion_d1}{should the d1 criterion be computed ? \emph{logical}.  Require \code{Wweight} to be computed.}  
  \item{criterion_entropy}{should the entropy criterion be computed ? \code{logical}.}  
  \item{criterion_Kalinsky}{should the Kalinsky criterion be computed ? \code{logical}.}  
  \item{criterion_Laboure}{should the Laboure criterion be computed ? \code{logical}.}
}
\references{
Chantal Revol-Muller, Francoise Peyrin, Yannick Carrillon and Christophe Odet. \emph{Automated 3D region growing algorithm based on an assessment function}. Pattern Recognition Letters, 23:137-150,2002.
}
\concept{calc.}
\keyword{functions,internal}
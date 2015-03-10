\name{calcHemi_cpp}
\title{Mid-saggital plan search}
\alias{calcHemi_cpp}
\description{
C++ function called by \code{\link{calcHemisphere}} to compute the objective function that evaluates the consistency of the mid-saggital plan. For internal use.
}
\usage{
calcHemi_cpp(px_hemiL, px_hemiR, sd_data, p, symetrie) 
}
\arguments{
  \item{px_hemiL}{the coordinates and the contrast of the observations in the left hemisphere for a given slice. \emph{matrix}.}
  \item{px_hemiR}{the coordinates and the contrast of the observations in the right hemisphere for a given slice. \emph{matrix}.}
  \item{sd_data}{the standard deviation of the parameter. \emph{numeric}.}
  \item{p}{the penalization factor. \emph{positive numeric}.}
  \item{symetrie}{the type of objective function. \code{TRUE} correspond to \code{"symmetry"} and \code{FALSE} to \code{"asymmetry"}.}
}
\concept{cpp.}
\keyword{functions,internal}
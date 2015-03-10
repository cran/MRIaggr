\name{calcRadius_cpp}
\title{Compute geometric caracteristics of a spatial group}
\alias{calcRadius_cpp}
\description{
C++ function called by \code{\link{calcGR}} to compute the barycenter of a spatial group. For internal use.
}
\usage{
calcRadius_cpp(coords, sample, threshold, subset_bary, trace)
}
\arguments{
  \item{coords}{the spatial coordinates of the observations. \emph{matrix} with a number of rows equal to the length of \code{sample}.}
  \item{sample}{the weight of each voxel in the computation of the barycenter. \emph{positive numeric}.}
  \item{threshold}{observations with a \code{sample} value below the value of \code{threshold} are discarded. \emph{numeric}.} 
  \item{subset_bary}{an indicator of the observations that should be kept ? \emph{logical vector}.}
  \item{trace}{should the radius of the spatial group be printed ? \emph{logical}.} 
}
\concept{cpp.}
\keyword{functions,internal}
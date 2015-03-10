\name{pointsHemisphere}
\title{Add the position of the mid-saggital plan}
\alias{pointsHemisphere}
\alias{pointsHemisphere,MRIaggr-method}
\description{
Add to an existing plot the mid-saggital plan. For internal use.
}
\usage{
\S4method{pointsHemisphere}{MRIaggr}(object,col="red",lwd=2,lty=1)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}.}
  \item{col}{the color with which the mid-saggital plan will be plotted. \emph{character}.}
  \item{lwd}{the line width. \emph{postive numeric}.}
  \item{lty}{the type of line used to represent the mid-saggital plan. \emph{numeric vector}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{lwd} and \code{lty} arguments can be found in \code{\link{par}}.
}
\concept{plot.}
\keyword{methods,internal}


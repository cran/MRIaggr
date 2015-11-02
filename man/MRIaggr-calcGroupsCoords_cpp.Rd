\name{calcGroupsCoords_cpp}
\title{Find spatial groups}
\alias{calcGroupsCoords_cpp}
\description{
C++ function called by \code{\link{calcGroupsCoords}} to compute the spatial groups. For internal use.
}
\usage{
calcGroupsCoords_cpp(coords_NNA, index_NNA, Neighborhood, coords_max, max_groups, verbose)
}
\arguments{
  \item{coords_NNA}{the spatial coordinates of the observations in C version (beginning at 0). \emph{matrix}.}
  \item{index_NNA}{the index of the coordinates in a array in C version (beginning at 0). \emph{numerical vector}.}
  \item{Neighborhood}{the type of neighbourhood. \emph{character}.}
  \item{coords_max}{the maximum coordinate in each dimension. \emph{numerical vector}}
  \item{max_groups}{the maximum number of groups. \emph{postive integer}.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\concept{cpp.}
\keyword{functions,internal}



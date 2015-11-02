\name{calcGroupsW_cpp}
\title{Find spatial groups}
\alias{calcGroupsW_cpp}
\description{
C++ function called by \code{\link{calcGroupsW}} to compute the spatial groups. For internal use.
}
\usage{
calcGroupsW_cpp(W, subset, max_groups)
}
\arguments{
  \item{W}{the neighbourhood matrix. \emph{dgCMatrix}.}
  \item{subset}{the subset of observations to use. \emph{positive integer vector}.}
  \item{max_groups}{the maximum number of groups. \emph{postive integer}.}
}
\concept{cpp.}
\keyword{functions,internal}



\name{filtrage2D_cpp}
\title{2D filtering}
\alias{filtrage2D_cpp}
\description{
C++ function call by \code{\link{calcFilter}} that performs two dimensional filtering. For internal use.
}
\usage{
filtrage2D_cpp(M_data, M_operateur, index_data, w_contrast, na_rm)
}
\arguments{
  \item{M_data}{matrix to which the filter will be applied.}
  \item{M_operateur}{the filter to be applied.}
  \item{index_data}{index of the non NA data.}
  \item{w_contrast}{should the influence of each neighbor be ponderated by the difference in signal with the considered observation ?}
  \item{na_rm}{should the observations with missing values in their neighborhood be removed ? Otherwise the ponderation is adjusted.}  
}
\concept{cpp.}
\keyword{functions,internal}
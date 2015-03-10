\name{selectVoxelSize}
\title{Extract the dimensions of a voxel}
\alias{selectVoxelSize}
\alias{selectVoxelSize,MRIaggr-method}
\description{
Extract the spatial dimensions of a voxel contained in a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectVoxelSize}{MRIaggr}(object,unit=TRUE)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{unit}{should the unit be returned ? \emph{logical}.}
}
\value{
  A \emph{data.frame} with one line and three or four columns containing the voxel size for each dimension and the size unit (if requested).
}
\examples{
data("MRIaggr.Pat1_red", package="MRIaggr")

## selection
selectVoxelSize(MRIaggr.Pat1_red)
}
\concept{select.}
\keyword{methods}


\name{selectVoxelDim}
\title{Extract the data dimension}
\alias{selectVoxelDim}
\alias{selectVoxelDim,Carto3D-method}
\alias{selectVoxelDim,MRIaggr-method}
\description{
Extract the spatial dimensions of the contrast data from a \code{\linkS4class{Carto3D}} or from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectVoxelDim}{Carto3D}(object)

\S4method{selectVoxelDim}{MRIaggr}(object)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.}
}
\value{
  A \emph{data.frame} with one line and three columns named \code{"i","j","k"} indicating the corresponding dimension.
}
\examples{
#### 1- Carto3D method ####
## load nifti files and convert them to Carto3D
path.nifti_files <- system.file("nifti",package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path.nifti_files,"TTP_t0"),format="nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")

## selection
selectVoxelDim(Carto3D.Pat1_TTP_t0)

#### 2- MRIaggr method ####
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## selection
selectVoxelDim(MRIaggr.Pat1_red)
}
\concept{select.}
\keyword{methods}


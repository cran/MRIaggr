\docType{class}
\name{Carto3D-class}
\title{class "Carto3D"}
\alias{Carto3D-class}
\alias{Carto3D}
\description{
  Patient-specific storage of univariate contrast data.
}
\arguments{
  \item{identifier}{the patient identifier. \emph{character}.}  
  \item{parameter}{the name of the contrast parameter. \emph{character}.}   
  \item{contrast}{the contrast and the spatial coordinates of each voxel. \emph{data.frame} with four columns named \code{"i"} \code{"j"} \code{"k"} and the name of the contrast parameter.}  
  \item{voxelDim}{the spatial dimensions of the lattice containing the observations. \emph{data.frame} with 1 row and 3 columns named \code{"i"} \code{"j"} \code{"k"}.}  
  \item{default_value}{the reference values of the contrast parameters (e.g. the background value). \emph{character}.}   
}

\section{S4 methods}{
\describe{
\item{\code{selectIdentifier}}{return the slot identifier.}
\item{\code{selectParameter}}{return the slot parameter.}
\item{\code{selectContrast}}{return the slot carto3D.}
\item{\code{selectVoxelDim}}{return the slot voxelDim.}
\item{\code{selectDefault_value}}{return the slot default_value.}

\item{\code{multiplot}}{display the contrast data.}

\item{\code{constMRIaggr}}{aggregate several Carto3D objects into a single objet MRIaggr.}

\item{\code{initNum}}{check the validity of the num argument.}
}
}
\seealso{
\code{\linkS4class{MRIaggr}} class that aggregates several \code{carto3D} objects.
\code{\link{Carto3D2MRIaggr}} to convert \code{carto3D} objects into \code{\linkS4class{MRIaggr}} objects.
}
\keyword{classes}

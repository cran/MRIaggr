\name{calcW}
\title{Compute the neighborhood matrix}
\alias{calcW}
\alias{calcW,data.frame-method}
\alias{calcW,MRIaggr-method}
\description{
Computate a neighborhood matrix using spatial coordinates.
}
\usage{
\S4method{calcW}{data.frame}(object,distband,method="euclidean",upper=NULL,
         format="dgCMatrix",row.norm=FALSE,spatial_res=rep(1,ncol(object)))

\S4method{calcW}{MRIaggr}(object,distband,spatial_res=c(1,1,1),num=NULL,
         hemisphere="both",subset=NULL,upper=TRUE,format="dgCMatrix",row.norm=FALSE,
         trace=TRUE,update.object=FALSE,overwrite=FALSE)
}
\arguments{
  \item{object}{a \emph{data.frame} containing the coordinates of the observations or an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{distband}{only distances smaller than delta are recorded. \code{postive numeric}. REQUIRED.}
  \item{method}{the distance measure to be used. \emph{character}. This must be one of \code{"euclidean"}, \code{"maximum"}, \code{"minkowski"} or \code{"greatcircle"}.}  
  \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to use. \emph{character}.}
  \item{subset}{the subset of observations to use. \emph{positive integer vector} or \code{NULL} leading to use all observations.}
  \item{upper}{should the entire matrix (\code{NULL}) or only the upper-triagonal (\code{TRUE}) or only the lower-triagonal (\code{FALSE}) values be calculated ?}
  \item{format}{the format of the output. Could be \code{"spam"} or \code{"dgCMatrix"}.}
  \item{row.norm}{should the resulting matrix be row-normalized ? \code{TRUE/FALSE}.}
   
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting neighborhood matrix be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if a neighborhood matrix is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{distband} argument can be found in \code{\link{nearest.dist}}.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.

The \code{row.norm} argument is ignored if \code{format} is set to \code{"spam"}.

FUNCTION: \cr
This function relies on the \code{nearest.dist} function of the \emph{spam} package.

Each of the \code{num}, \code{hemisphere} and \code{subset} argument define a subset of the total set of observations. 
It is the intersection of all these three subsets that is extracted.
}
\value{
  A \emph{spam} or \emph{dgCMatrix} object.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

#### 1- data.frame method ####
coords <- selectCoords(MRIaggr.Pat1_red,num=1:3,hemisphere="lesion")

## full W 
W <- calcW(object=coords,distband=sqrt(2))
W[1:10,1:10]
table(rowSums(W)    )

## full W normalized by row
W <- calcW(object=coords,distband=sqrt(2),row.norm=TRUE)
W[1:10,1:10]
table(rowSums(W)    )

## upper W 
W <- calcW(object=coords,distband=sqrt(2),upper=TRUE)
W[1:10,1:10]

#### 2- data.frame method ####

## compute W (regular lattice)
W <- calcW(MRIaggr.Pat1_red,distband=sqrt(2),upper=NULL,num=1:3,hemisphere="lesion")
table(rowSums(W>0))

## compute W (irregular lattice)
W <- calcW(MRIaggr.Pat1_red,distband=sqrt(2*1.875^2),upper=NULL,num=1:3,hemisphere="lesion",
           spatial_res=c(1.875,1.875,6))
table(rowSums(W>0))
}
\concept{calc.}
\keyword{methods}


                        
                        
                      
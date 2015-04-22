\name{affectDescStats}
\title{Affect non standard elements}
\alias{affectDescStats<-}
\alias{affectDescStats<-,MRIaggr-method}
\usage{
\S4method{affectDescStats}{MRIaggr}(object,name,overwrite=FALSE,trace=TRUE) <- value
}
\description{
Affect non standard elements to a \code{\linkS4class{MRIaggr}} object.
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{value}{any R object. REQUIRED.}
  \item{name}{the name of the element storing \code{value}. \emph{character}. REQUIRED.}
  \item{overwrite}{if an element with the same name is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.}  
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
FUNCTION: \cr
Contrary to all other \code{affect.} methods that impose restrictions on the objects that can be affected, this function enable to affect freely a R element.
However these elements will have no interaction with the methods of this package. 

There is an exception for the element \code{"W_euclidean"} that is a reserved space for a neighboring matrix. It must be a \code{dgCMatrix} with row size and column size corresponding to the number of observations.
}
\seealso{
 \code{\link{calcW}} to compute the neighboring matrix. \cr
 \code{\link{selectDescStats}} to extract non standard elements. 
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## affect a vector
affectDescStats(MRIaggr.Pat1_red,name="spatial_res") <- c(1.875,1.875,6)

## select the corresponding element
selectDescStats(MRIaggr.Pat1_red,"spatial_res")

## some calc methods automatically save results in the ls_descStats slot
# find spatial groups 
calcGroupsMask(MRIaggr.Pat1_red,mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
               W.range=6,W.spatial_res=selectDescStats(MRIaggr.Pat1_red,"spatial_res"),
               update.object=TRUE,overwrite=TRUE)

# extract spatial groups
selectDescStats(MRIaggr.Pat1_red,"GroupsLesion")

}
\concept{affect.}
\keyword{methods}

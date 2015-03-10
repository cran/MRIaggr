\name{selectDescStats}
\title{Extract non-stardard elements}
\alias{selectDescStats}
\alias{selectDescStats,MRIaggr-method}
\description{
Extract elements in the \code{ls_descStats} slot of a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectDescStats}{MRIaggr}(object,name=NULL,
         subset_W=NULL,hemisphere="both",num=NULL,slice_var="k")
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{name}{the name of the element to select. \emph{character} or \code{NULL} leading to select all available elements.}
  \item{subset_W}{the subset of observations to select. \emph{integer vector} or \code{NULL} leading to select all the observations.}
  \item{hemisphere}{the hemisphere to use. \emph{character}.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{slice_var}{the type of slice to extract. \code{"i"} for sagittal, \code{"j"} for coronal and \code{"k"} for transverse. \emph{character}.}  
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.

FUNCTION: \cr
Arguments \code{index.W}, \code{hemisphere}, \code{num}, \code{format.W}, \code{style.W} are active only if \code{name} equal \code{"W_euclidien"}. \cr
If \code{index.W} is not \code{NULL} then arguments \code{hemisphere} and \code{num} are ignored. \cr
Each of the \code{num}, \code{hemisphere} and \code{subset} argument define a subset of the total set of observations. 
It is the intersection of all these three subsets is extracted.
}
\seealso{ 
 \code{\link{affectDescStats<-}} to affect elements in the \code{ls_descStats} slot. \cr
 \code{\link{calcHemisphere}} to identify the hemispheres. \cr
 \code{\link{affectHemisphere<-}} and \code{\link{affectContrast<-}} to affect hemispheres obtained from an external source.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

calcGroupsMask(MRIaggr.Pat1_red,mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
               distband=6,
               spatial_res=c(1.875,1.875,6),
               update.object=TRUE,overwrite=TRUE)

## select all elements in the slot @ls_descStats
ls_descStats <- selectDescStats(MRIaggr.Pat1_red)
names(ls_descStats)

## get the name of all elements present in the slot @ls_descStats
selectParameter(MRIaggr.Pat1_red,type="ls_descStats")

## select a specific element 
res <- selectDescStats(MRIaggr.Pat1_red,name="GroupsLesion")

## compute and affect a neighborhood matrix
calcW(MRIaggr.Pat1_red,distband=3,update.object=TRUE,overwrite=TRUE,
      spatial_res=c(1.875,1.875,6))

## select the neighborhood matrix for a subset of observations
res <- selectDescStats(MRIaggr.Pat1_red,name="W_euclidean",hemisphere="lesion",num=3)
dim(res)
table(rowSums(res>0))

res <- selectDescStats(MRIaggr.Pat1_red,name="W_euclidean",subset_W = 1:10)
res
}
\concept{select.}
\keyword{methods}


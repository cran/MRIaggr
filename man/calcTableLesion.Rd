\name{calcTableLesion}
\title{Vertical distribution of the lesion}
\alias{calcTableLesion}
\alias{calcTableLesion,MRIaggr-method}
\description{
Compute the number of lesion observations at each slice.
}
\usage{
\S4method{calcTableLesion}{MRIaggr}(object,maskN,mask=NULL,as.logical=FALSE,
         trace=TRUE,update.object=FALSE,overwrite=FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{maskN}{the binary contrast parameter indicating the lesion. \emph{character}. REQUIRED.}
  \item{mask}{the binary contrast parameter indicating the brain. \emph{character} or \code{NULL} if no mask is available.}
  \item{as.logical}{should \code{mask} and \code{maskN} values be converted to logical ? \emph{logical}.}
 
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting lesion table be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if a lesion table is already stored in \code{object@table_lesion}, can it be overwritten ? \emph{logical}.}  
}
\value{
  A \emph{data.frame} containing the number of observation within the mask at each slice.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## compute table
res <- calcTableLesion(MRIaggr.Pat1_red,maskN=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
                       as.logical=TRUE,update.object=TRUE,overwrite=TRUE)

## extract table
res <- selectTable(MRIaggr.Pat1_red,"lesion")
}
\concept{calc.}
\keyword{methods}


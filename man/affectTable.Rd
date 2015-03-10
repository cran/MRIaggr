\name{affectTable}
\title{Affect volumic information}
\alias{affectTable<-}
\alias{affectTable<-,MRIaggr-method}
\description{
Affect volumic information to a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{affectTable}{MRIaggr}(object,type,overwrite=FALSE,trace=TRUE) <- value
}

\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{value}{the volumic information. \emph{data.frame}. REQUIRED.}
  \item{type}{the type of volumic information. Can be \code{"lesion"} \code{"reperfusion"} \code{"hypoperfusion"}. REQUIRED.}
  \item{overwrite}{if tables are already stored in \code{object@normalization}, can they be overwritten ? \emph{logical}.}  
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
ARGUMENTS : \cr
The validity of the \code{value} object is not checked. 
A valid format should match the result of the \code{calcTableHypoReperf} and \code{calcTableLesion} functions. 
}
\seealso{
 \code{\link{calcTableHypoReperf}} to compute the hypoperfusion and reperfusion tables. \cr
 \code{\link{calcTableLesion}} to compute the lesion table. \cr
 \code{\link{selectTable}} to extract the tables. \cr
 \code{\link{plotTableLesion}} to display the lesion volume by slice.
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

#### 1- lesion ####
## manual affectation
maskN <- c("MASK_DWI_t0","MASK_T2_FLAIR_t2")
resTable <- calcTableLesion(MRIaggr.Pat1_red,maskN=maskN,as.logical=TRUE)
affectTable(MRIaggr.Pat1_red,type="lesion",overwrite=TRUE) <- resTable

## automatic affectation
resTable <- calcTableLesion(MRIaggr.Pat1_red,maskN=maskN,
                            as.logical=TRUE,update.object=TRUE,overwrite=TRUE)

## display
selectTable(MRIaggr.Pat1_red,type="lesion")

#### 2- hypoperfusion and reperfusion ####
## manual affectation
resTable <- calcTableHypoReperf(MRIaggr.Pat1_red,param=c("TTP","MTT"),time=c("t0","t1"))
affectTable(MRIaggr.Pat1_red,type="hypoperfusion",overwrite=TRUE) <- resTable$volume_hypo
affectTable(MRIaggr.Pat1_red,type="reperfusion",overwrite=TRUE) <- resTable$volume_reperf

## automatic affectation
resTable <- calcTableHypoReperf(MRIaggr.Pat1_red,param=c("TTP","MTT"),time=c("t0","t1"),
                                update.object=TRUE,overwrite=TRUE)

## display
selectTable(MRIaggr.Pat1_red,type="reperfusion")
}
\concept{affect.}
\keyword{methods}


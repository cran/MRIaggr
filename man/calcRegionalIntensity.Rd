\name{calcRegionalContrast}
\title{Compute regional contrast parameters}
\alias{calcRegionalContrast}
\alias{calcRegionalContrast,MRIaggr-method}
\description{
Compute the regional contrast parameters.
}
\usage{
\S4method{calcRegionalContrast}{MRIaggr}(object,param,num=NULL,hemisphere="both",W="ifany",
         spatial_res=c(1,1,1),distband,distband_EDK,power_EDK=2,diagonal=FALSE,
         trace=TRUE,name_newparam=paste(param,"regional",sep="_"),
         update.object=FALSE,overwrite=FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter(s) from which the regional parameter(s) will be computed. \emph{character vector}. REQUIRED.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to use. \emph{character}.}
  \item{W}{the neighborhood matrix. \emph{dgCMatrix} or \code{"ifany"} leading to use the neighborhood matrix stored in the \code{object} if possible, else to compute this matrix.}
  \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.}
  \item{distband}{the neighborhood range. \emph{postive numeric}. Required only if \code{W} have to be computed.}
  \item{distband_EDK}{the distband of the kernel. \emph{postive numeric}. REQUIRED.}
  \item{power_EDK}{the power of the kernel. \emph{postive numeric}.}
  \item{diagonal}{should the diagonal be added to the neighborhood matrix ? \emph{logical}.}
  
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
  \item{name_newparam}{the name of the new parameters. \emph{character vector}.}
  \item{update.object}{should the resulting regional parameters be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if contrast parameters with the same names are already stored in \code{object} can they be overwritten ? \emph{logical}.}  
}
\details{
ARGUMENTS : \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere} arguments can be found in the details section of \code{\link{selectContrast}}.

Information about \code{distband_EDK} and \code{power_EDK} arguments can be found in the details section of \code{\link{EDK}}.
}
\value{
  A \emph{data.frame} containing in columns the regional parameters.
}
\seealso{
\code{\link{selectContrast}} to select the regional parameter(s). \cr
\code{\link{calcW}} to compute the neighboring matrix.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## compute regional values
\dontshow{
res  <- calcRegionalContrast(MRIaggr.Pat1_red,param=c("T2_FLAIR_t2","T1_t0"),
                        spatial_res=c(1.875,1.875,1.875),distband=6,distband_EDK=1.875)
}
\dontrun{
res  <- calcRegionalContrast(MRIaggr.Pat1_red,param=c("T2_FLAIR_t2","T1_t0"),
                        spatial_res=c(1.875,1.875,1.875),distband=6,distband_EDK=1.875,
                        update.object=TRUE,overwrite=TRUE)

## display
par(mfrow=c(2,4),mar=rep(1.5,4),mgp=c(2,0.5,0))
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",num=1:3,
             window=NULL,main="raw - slice ")
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2_regional",num=1:3,
             window=NULL,main="regional - slice ")
}
}
\concept{calc.}
\keyword{methods}


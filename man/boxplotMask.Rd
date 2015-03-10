\name{boxplotMask}
\title{Boxplot spatial group characteristics}
\alias{boxplotMask}
\alias{boxplotMask,MRIaggr-method}
\description{
Display a boxplot of the contrast parameter values inside and outside a spatial group.
}
\usage{
\S4method{boxplotMask}{MRIaggr}(object,param,mask,num=NULL,hemisphere="both",
         norm_mu=FALSE,norm_sigma=FALSE,scale=TRUE,as.logical=FALSE,
         window=FALSE,ylim=NULL,col=c("white","purple"),main=NULL,mgp=c(2,0.5,0),
         x.legend="topright",y.legend=NULL,cex.legend=0.8,
         filename=paste(object@identifier,"boxplotMask",sep="_"),
         width=1000,height=700,path=NULL,unit="px",res=NA)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
 
  \item{param}{the contrast parameter(s) associated with the lesion mask(s). \emph{character vector}. REQUIRED.}
  \item{mask}{the binary contrast parameter(s) defining the spatial group(s). \emph{character vector}. REQUIRED.}
  \item{num}{the slices to consider. \emph{numeric vector} or \code{NULL}. REQUIRED.}
  \item{hemisphere}{the hemisphere to consider. \emph{character}.}
  \item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}.}
  \item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}.}

  \item{scale}{should the contrast parameters be scaled ? \emph{logical}.}
  \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.}
   
  \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.}
  \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.}
  \item{col}{the colors of the boxplots for observations inside and outside the mask(s). \emph{character vector of size 2}.}
  \item{main}{an overall title for the plot. \emph{character}.}
  \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.}
  \item{x.legend}{the x coordinates of the legend. \emph{numeric} or \emph{character}.}
  \item{y.legend}{the y coordinates of the legend. \emph{numeric} or \emph{character}.}
  \item{cex.legend}{the expansion factor of the legend. \emph{positive numeric}.}		
 
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{width}{the width of the device used to export the plot. \emph{postive numeric}.}
  \item{height}{the height of the device used to export the plot. \emph{postive numeric}.}
  \item{path}{the directory where the plot file will be created. \emph{character}.}
  \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.}
  \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer} or \code{NA}.}
}
\value{
  None.
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere}, \code{norm_mu} and \code{norm_sigma} arguments can be found in the details section of \code{\link{selectContrast}}.

Information about the \code{window}, \code{filename}, \code{width}, \code{height}, \code{path}, \code{unit} and \code{res} arguments can be found in the details section of \code{\link{initWindow}}.

Information about the \code{mgp} argument can be found in \code{\link{par}}.

Information about the \code{x.legend}, \code{y.legend}, \code{cex.legend} arguments can be found in \code{\link{legend}} (\code{cex.legend} is the \code{cex} argument of \code{legend}).
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## display
boxplotMask(MRIaggr.Pat1_red,param=c("DWI_t0","TTP_t0","MTT_t0"),mask="MASK_T2_FLAIR_t2",
             as.logical=TRUE)
}
\concept{plot.}
\keyword{methods}


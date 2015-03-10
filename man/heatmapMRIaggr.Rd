\name{heatmapMRIaggr}
\title{Correlation between contrast parameters}
\alias{heatmapMRIaggr}
\alias{heatmapMRIaggr,MRIaggr-method}
\description{
Display a correlation map of the contrast parameters.
}
\usage{
\S4method{heatmapMRIaggr}{MRIaggr}(object,param,num=NULL,hemisphere="both",
         scale=TRUE,method="pearson",points.values=TRUE,type="image",digit=3,
         breaks=NULL,window=FALSE,col=cm.colors(256),main=NULL,
         mgp=c(2,0.5,0),mar=c(4,4,1,6),las=1,cex.axis=1,
         filename=paste(object@identifier,"heatmapMRIaggr",sep="_"),
         width=1000,height=700,path=NULL,unit="px",res=NA)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameters used to compute the correlations. \emph{character vector}. REQUIRED.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to use. \emph{character}. See the details section of \code{\link{selectContrast}}.}
  \item{scale}{should the contrast parameters be scaled ? \emph{logical}.}

  \item{method}{the correlation coefficient which is to be computed. Can be \code{"pearson"}, \code{"kendall"} or \code{"spearman"}.}
  \item{points.values}{should the correlation values be printed on the plot ? \emph{logical}.}
  \item{type}{the type of plot to display. Any of \code{"image"} or \code{"image.plot"} or \code{FALSE} meaning no plot.}
  \item{digit}{the number of decimal places to use for the labels. \emph{positive integer}.}
  \item{breaks}{the break points to use to generate the color intervals. \emph{numeric vector} or \code{NULL} leading to automatic breakpoints generation.}

  \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.}
  \item{col}{the colors with which the correlations will be displayed. \emph{character vector}.}
  \item{main}{an overall title for the plot. \emph{character}.}
  \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.}
  \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.}
  \item{las}{the style of the axis labels. Any of \code{0}, \code{1}, \code{2} or \code{3}.}
  \item{cex.axis}{the magnification to be used for axis annotation relative to the current setting of cex. \emph{positive numeric}.}

  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{width}{the width of the device used to export the plot. \emph{postive numeric}.}
  \item{height}{the height of the device used to export the plot. \emph{postive numeric}.}
  \item{path}{the directory where the plot file will be created. \emph{character}.}
  \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.}
  \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.

If \code{breaks} is not \code{NULL}, it must be of length \code{length(col)+1}.

Information about the \code{window}, \code{filename}, \code{width}, \code{height}, \code{path}, \code{unit} and \code{res} arguments can be found in the details section of \code{\link{initWindow}}.

Information about the \code{mar}, \code{las} and \code{mgp} arguments can be found in \code{\link{par}}.
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red",package="MRIaggr")

## pearson
heatmapMRIaggr(MRIaggr.Pat1_red,param=c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las=1,type="image",cex=0.75,
           breaks=seq(-1,1,length.out=51),
           col=cm.colors(50))  	
## spearman
heatmapMRIaggr(MRIaggr.Pat1_red,param=c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las=1,type="image",cex=0.75,method="spearman",
           breaks=seq(-1,1,length.out=51),
           col=cm.colors(50))  


## spearman with legend
  heatmapMRIaggr(MRIaggr.Pat1_red,param=c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las=1,type="image.plot",cex=0.75,method="spearman",
           breaks=seq(-1,1,length.out=51),
           col=cm.colors(50))  

}
\concept{plot.}
\keyword{methods}


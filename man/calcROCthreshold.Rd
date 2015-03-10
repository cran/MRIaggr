\name{calcROCthreshold}
\title{ROC analysis}
\alias{calcROCthreshold}
\alias{calcROCthreshold,MRIaggr-method}
\description{
Perform a ROC analysis on a continuous variable for identifying a binary outcome.
}
\usage{
\S4method{calcROCthreshold}{MRIaggr}(object,param,mask,as.logical=FALSE,digit=2,
         plot="ROC_Youden",digit.plot=3,
         window=FALSE,filename=paste(object@identifier,"calcROCthreshold",plot,sep="_"),
         width=1000,height=700,path=NULL,unit="px",res=NA,
         trace=TRUE,update.object=FALSE,overwrite=FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter(s) that should be used to identify the observations inside the mask. \emph{character vector}. REQUIRED.}
  \item{mask}{the binary contrast parameter that will be used as the outcome in the ROC analysis. \emph{character vector}. REQUIRED.}
  \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.}
  \item{digit}{the number of decimal places to use when generating the thresholds. \emph{positive integer}.}
  \item{plot}{the type of the graphic to display? \emph{character} or \code{FALSE}. See the details section.}
  \item{digit.plot}{the number of decimal places to use for the legend. \emph{positive integer}.}
 
  \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.}
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{width}{the width of the device used to export the plot. \emph{postive numeric}.}
  \item{height}{the height of the device used to export the plot. \emph{postive numeric}.}
  \item{path}{the directory where the plot file will be created. \emph{character}.}
  \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.}
  \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.}
 
  \item{update.object}{should the resulting threshold analysis be stored in \code{object@ls_descStats} ? \emph{logical}.}
  \item{overwrite}{if a threshold analysis is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.}  
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{window}, \code{filename}, \code{width}, \code{height}, \code{path}, \code{unit} and \code{res} arguments can be found in the details section of \code{\link{initWindow}}.

Possible values for \code{plot} are:
\itemize{
\item \code{"ROC_Youden"} : display the ROC curve with the optimal threshold according to the Youden index.
\item \code{"ROC_prev"} : display the ROC curve with the optimal threshold according to the utility function.
\item \code{"boxplot_Youden"} : display a boxplot of the contrast parameter for each outcome with the optimal threshold according to the Youden index.
\item \code{"boxplot_prev"} : display a boxplot of the contrast parameter for each outcome with the optimal threshold according to the utility function.
\item \code{FALSE} : no graphic is displayed.
}

}
\value{
   An \emph{data.frame} containing for each mask the AUC and AUPRC value, the optimal threshold and the corresponding sensitivity and specificity for the Youden criteria and a utility function taking into account the prevalence.
}
\seealso{
\code{\link{selectDescStats}} to select the mask caracteristics.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")
  
## ROC analysis
res <- calcROCthreshold(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),
                        mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),as.logical=TRUE)
  
res <- calcROCthreshold(MRIaggr.Pat1_red,param=c("DWI_t0","T2_FLAIR_t2"),
                        mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),as.logical=TRUE,
                        plot="boxplot_Youden",
						update.object=TRUE,overwrite=TRUE)
						
selectDescStats(MRIaggr.Pat1_red,"Mask_threshold")
}
\concept{calc.}
\keyword{methods}


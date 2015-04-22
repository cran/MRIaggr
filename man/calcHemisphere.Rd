\name{calcHemisphere}
\title{Find the mid-saggital plan}
\alias{calcHemisphere}
\alias{calcHemisphere,MRIaggr-method}
\description{
Find a plane that distinguish the two cerebral hemispheres.
}
\usage{
\S4method{calcHemisphere}{MRIaggr}(object,param,num=NULL,p=1,subset=NULL,
         penalty="symmetry",mask=NULL,as.logical=FALSE,i_test=5,j_test=5,
         angle_test=5,unit_angle="radian",n.points=100,plot=TRUE,
         window=FALSE,filename=paste(object@identifier,"_calcHemisphere",sep=""),
         width=1000,height=700,path=NULL,unit="px",res=NA,
         trace=TRUE,update.object=FALSE,overwrite=FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter that should be used to distinguish the two hemispheres. \emph{character}. REQUIRED.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{p}{the type of distance for the penalization. \emph{positive numeric}.}
  \item{subset}{the subset of observations to use. \emph{positive integer vector} or \code{NULL} leading to use all observations.}
  \item{penalty}{the type of objective function. Can be \code{"symmetry"} or \code{"asymmetry"}.}
  \item{mask}{the binary contrast parameter(s) indicating the lesion. \emph{character vector} or \code{NULL} if no mask is available.}
  \item{as.logical}{should \code{mask} be convert to logical ? \emph{logical}.}

  \item{i_test}{the abscissa or the number of abscissa to test. \emph{numeric vector} or \emph{positive integer}.}
  \item{j_test}{the ordinates or the number of ordinates to test. \emph{numeric vector} or \emph{positive integer}.}
  \item{angle_test}{the angle or the number of angle to test. \emph{numeric vector} or \emph{postive integer}.}
  \item{unit_angle}{the unit in which the angle is given. Can be \code{"radian"} or \code{"degree"}.}
  \item{n.points}{the number of points that represent the mid-saggital plan to computed. \code{positive integer}.}

  \item{plot}{should the results be plotted ? \emph{logical}.}
  \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.}
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{width}{the width of the device used to export the plot. \emph{postive numeric}.}
  \item{height}{the height of the device used to export the plot. \emph{postive numeric}.}
  \item{path}{the directory where the plot file will be created. \emph{character}.}
  \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.}
  \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.}
  
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting midplane be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if a midplane is already stored in \code{object@midplane}, can it be overwritten ? \emph{logical}.}  
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Setting \code{p} to \code{1} leads to use the absolute difference ; setting \code{p} to \code{2} leads to use the euclidean distance.

Information about the \code{window}, \code{filename}, \code{width}, \code{height}, \code{path}, \code{unit} and \code{res} arguments can be found in the details section of \code{\link{initWindow}}.

FUNCTION: \cr
This function seeks the plane that minimize the difference between controlateral values of the two hemispheres or maximize the similarity between the two hemispheres.
There are 3 degree of freedom : two for the position of the center and one for the angle.
The separation between the hemisphere is assumed to be identical for all slices. 
From our experience, using an objective function based on symmetry gives better results compared to asymmetry.
}
\seealso{
\code{\link{selectParameter}} to select the midplane.
}
\value{
 An \emph{list} containing :
  \itemize{
  \item\code{[[penalty]]} : an \emph{array} containing the objective function for the various configurations.
  \item\code{[[nb]]} : an \emph{array} containing the number of observations used to compute the penalty function for the various configurations.
  \item\code{[[moy]]} : an \emph{array} containing the mean value of the objective function for the various configurations.
  \item\code{[[optimum]]} : the parameters of the optimal midplane. \code{numeric vector}
  \item\code{[[midplane]]} : the position of the midplane points. \emph{matrix}.
  \item\code{[[data]]} : the position of the observations relatively to the mid-saggital plane. \emph{matrix}.
 \item\code{[[cv]]} : Was the optimum reached inside the parameter space and not at a border ? \emph{logical}.
  }
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

\dontrun{
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)
}
\dontshow{
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",num=1,
                      mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),i_test=1,j_test=1,angle_test=1,
                      trace=TRUE,update.object=TRUE,overwrite=TRUE)
}

## display the mid-saggital plan
multiplot(MRIaggr.Pat1_red,param="T2_GRE_t0",num=3,legend=FALSE,
             midplane=TRUE,main="original coordinates - slice ")

## display
multiplot(selectContrast(MRIaggr.Pat1_red,param=c("i_hemisphere","j_hemisphere","k")),
             contrast=selectContrast(MRIaggr.Pat1_red,param="T2_GRE_t0"),num=3,
             index1=cbind(0,seq(-50,50),3),main="new coordinates - slice ",legend=FALSE)


## compute the mid-saggital plan and mark lesion
\dontrun{
res <- calcHemisphere(MRIaggr.Pat1_red,param="T2_GRE_t0",
         mask=c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
         trace=TRUE,update.object=TRUE,overwrite=TRUE)
}
}
\concept{calc.}
\keyword{methods}


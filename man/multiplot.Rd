\name{multiplot}
\title{Slice by slice display}
\alias{multiplot}
\alias{multiplot,data.frame-method}
\alias{multiplot,Carto3D-method}
\alias{multiplot,MRIaggr-method}
\description{
Make a slice by slice display of a \emph{data.frame}, a \code{\linkS4class{Carto3D}} or a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{multiplot}{data.frame}(object,contrast=NULL,num=NULL,
         index1=NULL,index2=NULL,index3=NULL,breaks=50,type.breaks="range",
         palette="terrain.colors",col=NULL,pch=NULL,cex=1,col.NA="lightyellow",pch.NA=8,
         xlim=NULL,ylim=NULL,axes=TRUE,window=FALSE,legend=TRUE,mfrow=NULL,
         mar=rep(1.5,4),mgp=c(2,0.5,0),pty=NULL,asp=1,bg="lightblue",xlab="",ylab="",
         main="slice ",num.main=TRUE,cex.main=1.5,quantiles.legend=TRUE,digit.legend=3,
         cex.legend=1.5,mar.legend=c(2,7,2,2),main.legend=names(contrast),
         filename="multiplot",width=1000,height=500,path=NULL,unit="px",res=NA)
		 
\S4method{multiplot}{Carto3D}(object,num=NULL,
         breaks=50,type.breaks="range",palette="terrain.colors",col=NULL,pch=NULL,cex=1,
         col.NA="lightyellow",pch.NA=8,xlim=NULL,ylim=NULL,axes=TRUE,window=FALSE,
         legend=TRUE,mfrow=NULL,mar=rep(1.5,4),mgp=c(2,0.5,0),pty=NULL,asp=1,
         bg="lightblue",xlab="", ylab="",main=NULL,num.main=TRUE,cex.main=1.5,
         quantiles.legend=TRUE,digit.legend=3,cex.legend=1.5,mar.legend=c(2,7,2,2),
         filename="multiplot",width=1000,height=500,path=NULL,unit="px",res=NA)

\S4method{multiplot}{MRIaggr}(object,param,num=NULL,
         index1=NULL,index2=NULL,index3=NULL,midplane=FALSE,slice_var="k",
         hemisphere="both",norm_mu=FALSE,norm_sigma=FALSE,as.logical=FALSE,
         breaks=50,type.breaks="range",palette="terrain.colors",col=NULL,pch=NULL,cex=1,
         col.NA="lightyellow",pch.NA=8,col.midplane="red",xlim=NULL,ylim=NULL,axes=TRUE,
         window=FALSE,legend=TRUE,mfrow=NULL,mar=rep(1.5,4),mgp=c(2,0.5,0),pty=NULL,
         asp=1,bg="lightblue",xlab="", ylab="",main=NULL,num.main=TRUE,cex.main=1.5,
         quantiles.legend=TRUE,digit.legend=3,cex.legend=1.5,
         mar.legend=c(2,7,2,2),main.legend = param,
         filename="multiplot",width=1000,height=500,path=NULL,unit="px",res=NA)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}} or \code{\linkS4class{Carto3D}} or a 3 column \code{data.frame} containing the coordinates of the observations in columns. REQUIRED.}
  \item{param}{the contrast parameter to display. \emph{character}. REQUIRED.}
  \item{contrast}{the intensities to display. \emph{numerical vector} or \code{NULL} leading to use the same color for all observations.}
  \item{num}{the slices to display. \emph{numeric vector} or \code{NULL}.}
  \item{index1,index2,index3}{the coordinates of additionnal points to display. \emph{data.frame} or \emph{list} or \code{NULL}.}
  \item{midplane}{should the mid-sagittal plan be displayed ? \emph{logical}.}
 
  \item{slice_var}{the type of view to use. \code{"i"} for sagittal view, \code{"j"} for coronal view and \code{"k"} for transverse view. \emph{character}.}  
  \item{hemisphere}{the hemisphere to display. \emph{character}.}
  \item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}.}
  \item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}.}
  \item{as.logical}{if a parameter is specified for the index arguments, should it be converted to logical ? \emph{logical}.}
   
  \item{breaks}{the break points or the number of breakpoints to use to generate the color intervals . \emph{numeric vector} or \emph{postive integer}.}
  \item{type.breaks}{should the break points be equally space according the range of data values (\code{"range"}), centered (\code{"range_center"}) or correspond to the quantile values (\code{"quantile"}) ?.}
  \item{palette}{the colors or the palette to use when associating colors to observation intensities. \emph{character vector} or \emph{character}.}
  \item{col}{the color with which the observations will be displayed. \emph{character vector} or \code{NULL} leading to determine the colors using the \code{palette} and \code{breaks} arguments.}
  \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer} or \code{NULL} leading to use the \code{image} function instead of \code{plot}.}
  \item{cex}{the expansion factor used to plot the observations. \emph{positive numeric}.}

  \item{col.NA}{the color to use to plot the NAs. \emph{character}.}
  \item{pch.NA}{the label to use to plot the NAs. \emph{postive integer}.}
  \item{col.midplane}{the color in which the mid-sagittal plan should appear. \emph{character}.}
  \item{xlim}{the x limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the x limits.}
  \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.}
  \item{axes}{should the axes be plotted ? \emph{logical}.}
 
  \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.}
  \item{legend}{how should the legend be displayed ? \emph{logical} or \code{NULL}.}
  \item{mfrow}{the division of the device in plot region. \emph{numeric vector of size 2} or \code{NULL} leading automatic adjustment.}
  \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.}
  \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.}
  \item{pty}{the type of plot region to be used. Can be \code{"s"} or \code{"m"}.}  
  \item{asp}{the aspect ratio y/x. \emph{numeric}.}
  \item{bg}{the color used for the background. \emph{character}.}

  \item{xlab}{a title for the x axis. \emph{character}.}
  \item{ylab}{a title for the y axis. \emph{character}.}
  \item{main}{an overall title for the plot. \emph{character}.}
  \item{num.main}{should the slice number be written over each plot. \emph{logical}.}
  \item{cex.main}{the expansion factor for the main title. \emph{numeric}.}

  \item{quantiles.legend}{should the quantiles of the data be displayed on the legend ? \emph{logical}.}
  \item{digit.legend}{the number of decimal places to use for the legend labels. \emph{integer}.}
  \item{cex.legend}{the expansion factor of the legend. \emph{positive numeric}.}
  \item{mar.legend}{the number of margin lines to be specified on the four sides of the legend. \emph{numeric vector  of size 4}.}
  \item{main.legend}{a main title for the legend. \emph{character}.}
 
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

Information about the \code{index1}, \code{index2} and \code{index3} arguments can be found in the details section of \code{\link{initIndex}} (argument \code{index}).

Information about the \code{hemisphere}, \code{norm_mu} and \code{norm_mu} arguments can be found in the details section of \code{\link{selectContrast}}.

Information about the \code{window}, \code{filename}, \code{width}, \code{height}, \code{path}, \code{unit} and \code{res} arguments can be found in the details section of \code{\link{initWindow}}.

Information about the \code{mfrow}, \code{mar}, \code{mar.legend}, \code{mgp} \code{pty}, \code{asp} and \code{bg} arguments can be found in \code{\link{par}}.

Information about the \code{palette} and \code{legend} arguments can be found in the details section of \code{\link{initCol}}.
}
\value{
  None.
}
\seealso{ 
 \code{\link{plotLesion3D}} for a 3D plot of the lesion. \cr
 \code{slices3d} of the \emph{misc3d} package for a more interactive 3D plot.
}
\examples{
#### 1- data.frame ####
## simulate 
n <- 10
Y <- rnorm(n^2)

## display
multiplot(object=data.frame(expand.grid(1:n,1:n),1),
             contrast=Y,window=FALSE)

## load a MRIaggr object
data(MRIaggr.Pat1_red, package="MRIaggr")

## select data
data <- selectContrast(MRIaggr.Pat1_red,param=c("DWI_t0","TTP_t0","MTT_t0","MASK_T2_FLAIR_t2"),
                    hemisphere="lesion",coords=TRUE)
					
## fit model
glm.1 <- glm(MASK_T2_FLAIR_t2 ~ DWI_t0 + TTP_t0 + MTT_t0,data=data,family=binomial(link="logit"))

## display fitted values
multiplot(object=data[,c("i","j","k")],
             contrast=predict(glm.1,type="response"),window=FALSE)

## display residuals
multiplot(object=data[,c("i","j","k")],num=3,
             contrast=predict(glm.1,type="response"),window=FALSE,
             index1=list(coords=data[data$MASK_T2_FLAIR_t2,c("i","j","k")],outline=TRUE)
             )
			  
#### 2- carto3D ####
## load NIFTI files and convert them to carto3D
path.Pat1 <- system.file("nifti",package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path.Pat1,"TTP_t0"),format="nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0,identifier="Pat1",param="TTP_t0")

## display
multiplot(Carto3D.Pat1_TTP_t0)
\dontrun{
multiplot(Carto3D.Pat1_TTP_t0,num=1:2)
multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=FALSE)
multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=FALSE,legend=FALSE)
multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=FALSE,legend=FALSE,main="",num.main=FALSE)
multiplot(Carto3D.Pat1_TTP_t0,num=1:2,axes=FALSE,main="",num.main=FALSE,
			 palette="gray.colors",breaks=seq(0,100))
}

#### 3- MRIaggr ####
## load a MRIaggr object
data(MRIaggr.Pat1_red, package="MRIaggr")

# display 3 slices 
multiplot(MRIaggr.Pat1_red,param="DWI_t0",              
             num=1:3)

\dontrun{
# display 3 slices with no axes and white background
multiplot(MRIaggr.Pat1_red,param="DWI_t0",
             num=1:3,axes=FALSE,bg="white")

# remove the legend
multiplot(MRIaggr.Pat1_red,param="DWI_t0",              
             num=1:3,legend=FALSE)
}

## display an set of points
# using a binary parameter stored in the object
multiplot(MRIaggr.Pat1_red,param="DWI_t0",              
             num=1:3,index1=list(coords="MASK_DWI_t0")
)

\dontrun{
# customize the display of the points
multiplot(MRIaggr.Pat1_red,param="DWI_t0",              
             num=1:3,index1=list(coords="MASK_DWI_t0",col="pink",pch=14)
)

# display only the edges of the set
multiplot(MRIaggr.Pat1_red,param="DWI_t0",num=3, legend=FALSE,  
             index1=list(coords="MASK_DWI_t0",outline=TRUE)
)
             

# specify the index of points using coordinates
coordsIndex <- data.frame(i=c(40,60),j=c(80,100),k=c(3,3))

multiplot(MRIaggr.Pat1_red,param="DWI_t0",num=3,legend=FALSE,       
             index2=list(coords=coordsIndex,col="black",pch=15,cex=4)
             )


# various possibilities for the display
multiplot(MRIaggr.Pat1_red,num=1:3,param="DWI_t0",
             legend=FALSE,window=FALSE)
multiplot(MRIaggr.Pat1_red,num=1:3,param="DWI_t0",
             legend=TRUE,window=FALSE)
multiplot(MRIaggr.Pat1_red,num=1:3,param="DWI_t0",
             legend=NULL,window=FALSE)
			 }
}
\concept{plot.}
\keyword{methods}


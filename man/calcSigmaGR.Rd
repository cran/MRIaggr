\name{calcSigmaGR}
\title{Automatic Growing Region algorithm}
\alias{calcSigmaGR}
\description{
Evaluate the quality of the Growing Region partition for several homogeneity parameters.
}
\usage{
calcSigmaGR (contrast,W,seed,sigma,range=c(-Inf,+Inf),range.seed=c(-Inf,+Inf),
         breaks=100,scale=FALSE,iter_max=100,sd.robust=FALSE,criterion_d1=FALSE,
         criterion_entropy=TRUE,criterion_Kalinsky=TRUE,criterion_Laboure=TRUE,
         trace=TRUE,mar=rep(2,4),mgp=c(2,0.75,0),main="",window=FALSE,
         filename="calcSigmaGR",width=1000,height=700,path=NULL,unit="px",res=NA)
}
\arguments{
  \item{contrast}{the contrast value of each observation. \emph{numeric vector}. REQUIRED.}
  \item{W}{the neighborhood matrix. \code{dgCMatrix}. REQUIRED.}
  \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \code{positive integer vector} or \code{logical vector}. REQUIRED.}
  \item{sigma}{the sequence of maximum admissible values for the group variability \code{positive numeric vector}. REQUIRED.}
  \item{range}{the range of acceptable contrast values for the growing region group. \code{numeric vector of size 2}.}
  \item{range.seed}{the range of acceptable contrast values for the seeds. \code{numeric vector of size 2}.}
 
  \item{breaks}{the break points or the number of break points to use to categorize the contrast distribution. \code{numeric vector} or \code{postive integer}.}
  \item{scale}{should the contrast be scaled ? \code{logical}.}  
  \item{iter_max}{the maximum number of iterations for the expansion of the growing region. \code{postive integer}.}  
  \item{sd.robust}{should the median absolute deviation be used to estimte the variability of the group contrast (\code{TRUE}), or the standard deviation (\code{FALSE}) ? \code{logical}.}
  
  \item{criterion_d1}{should the d1 criterion be computed ? \emph{logical}.  Require \code{Wweight} to be computed.}  
  \item{criterion_entropy}{should the entropy criterion be computed ? \code{logical}.}  
  \item{criterion_Kalinsky}{should the Kalinsky criterion be computed ? \code{logical}.}  
  \item{criterion_Laboure}{should the Laboure criterion be computed ? \code{logical}.}
  
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}  
  \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.}
  \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.}
  \item{main}{an overall title for the plot. \emph{character}.}

  \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \code{character}.}
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{width}{the width of the device used to export the plot. \emph{postive numeric}.}
  \item{height}{the height of the device used to export the plot. \emph{postive numeric}.}
  \item{path}{the directory where the plot file will be created. \emph{character}.}
  \item{unit}{the units in which \code{height} and \code{width} are given. \emph{character}.}
  \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{window}, \code{filename}, \code{width}, \code{height}, \code{path}, \code{unit} and \code{res} arguments can be found in the details section of \code{\link{initWindow}}.

Information about the \code{mar} and \code{mgp} arguments can be found in \code{\link{par}}.

FUNCTION: \cr
The Growing Region step has been proposed by (Revol et al. 1997) and (Revol et al. 2002).
}
\value{
  An \emph{list} containing :
  \itemize{
  \item\code{[[df.criterion]]} : the value of the clustering criteria (in columns) for each sigma value (in rows). \emph{numeric matrix}.
  \item\code{[[list.GR]]} : the list of the optimal GR sets, one for each clustering criteria.
  \item\code{[[best]]} : the optimal value of each clustering criteria. \emph{data.frame}.
 }
}
\references{
Chantal Revol and Michel Jourlin. \emph{A new minimum varance region growing algorithm for image segmentation}. Pattern Recognition Letters, 18(3):249-258,1997. \cr
Chantal Revol-Muller, Francoise Peyrin, Yannick Carrillon and Christophe Odet. \emph{Automated 3D region growing algorithm based on an assessment function}. Pattern Recognition Letters, 23:137-150,2002.
}
\examples{
\dontrun{
## load a MRIaggr object
data(MRIaggr.Pat1_red, package="MRIaggr")

## select data
carto <- selectContrast(MRIaggr.Pat1_red,num=3,param=c("TTP_t0","MASK_DWI_t0"),hemisphere="lesion")
coords <- selectCoords(MRIaggr.Pat1_red,num=3,hemisphere="lesion")
W <- calcW(coords,distband=sqrt(2))

indexN <- which(carto$MASK_DWI_t0==1)
seed <- indexN[which.max(carto[indexN,"TTP_t0"])]

## find optimal sigma
resGR_sigma <- calcSigmaGR(contrast=carto$TTP_t0, W=W,seed=seed,
                           sigma=seq(0.8,2,0.1),iter_max=50)
}
}
\concept{calc.}
\keyword{functions}
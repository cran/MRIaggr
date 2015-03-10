\name{calcAUPRC}
\title{Area under the PR curve}
\alias{calcAUPRC}
\description{
Compute the area under the precision recall curve by numerical integration.
}
\usage{
calcAUPRC(x,y,subdivisions=10000,performance=NULL)
}
\arguments{
  \item{x}{the biomarker values. \emph{numeric vector}. REQUIRED.}
  \item{y}{the class labels. \emph{numeric vector}, \emph{character vector} or \emph{logical vector}. REQUIRED.}
  \item{subdivisions}{the maximum number of subintervals used for the integration. \emph{positive integer}.}
  \item{performance}{an object of class \code{performance} can be supplied instead of arguments \code{x} and \code{y}.}  
}
\details{
ARGUMENTS: \cr
\code{y} must have exactly two levels.

If \code{performance} is set to \code{NULL}, the code{x} and \code{y} will be used to form the \code{performance} object.
}
\value{
  A \emph{numeric} between 0 and 1.
}
\examples{
#### 1- with MRIaggr data ####
## load a MRIaggr object
data(MRIaggr.Pat1_red,package="MRIaggr")

## select parameter and binary outcome
cartoT2 <- selectContrast(MRIaggr.Pat1_red,param="T2_FLAIR_t2")
cartoMASK <- selectContrast(MRIaggr.Pat1_red,param="MASK_T2_FLAIR_t2")

## compute AUPRC
T2.AUPRC <- calcAUPRC(x=cartoT2,y=cartoMASK)

## compute AUC
# if(require(pROC)){
# T2.AUC <- auc(roc(cartoMASK ~ cartoT2))
# }

## display
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",num=1,
             index1=list(coords="MASK_T2_FLAIR_t2",outline=TRUE)
             )

#### 2- with simulated data ####
n0 <- 1000
n1 <- c(10,100,1000)
for(iter_n in 1:length(n1)){
  X <- c(rnorm(n0,0),rnorm(n1[iter_n],2))
  Y <- c(rep(0,n0),rep(1,n1[iter_n]))
  print(calcAUPRC(X,Y))
}

## alternative way using a performance object
perfXY <- ROCR::performance(ROCR::prediction(X,Y),x.measure="rec",measure="prec")
calcAUPRC(performance=perfXY,subdivisions=10000)
}
\concept{const.}
\keyword{functions}


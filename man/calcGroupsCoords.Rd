\name{calcGroupsCoords}
\title{Compute spatial groups}
\alias{calcGroupsCoords}
\description{
Compute the spatial groups using the coordinates of the observations.
}
\usage{
calcGroupsCoords(coords,array=NULL,Neighborhood,max_groups=10000,trace=TRUE)
}
\arguments{
  \item{coords}{the spatial coordinates of the observations. \emph{data.frame}. REQUIRED.}
  \item{array}{alternative specification of the spatial coordinates using an array where the non-NA values indicates the points of interest. \emph{array} or \code{NULL} leading to consider the \code{coords} argument.}
  \item{Neighborhood}{the type of neighborhood. \emph{character}.}
  \item{max_groups}{the maximum number of groups. \emph{postive integer}.}
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
the \code{Neighborhood} argument can be a \emph{matrix} or an \emph{array} defining directly the neighborhood to use (i.e the weight of each neighbor) 
or a name indicating which type of neighborhood should be used (see the details section of \code{\link{initNeighborhood}}).
}
\value{
 An \emph{list} containing :
  \itemize{
  \item\code{[[ls.group]]} : \code{lists} of the observations of each spatial group.
  \item\code{[[df.group]]} : a \emph{data.frame} indicating the position and the group of each observation.
  \item\code{[[group_size]]} : the size of each spatial group. \emph{integer vector}.
  }
}
\seealso{
\code{\link{initFilter}} for various pre-stored filters.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select data
MASK_DWI_t0 <- selectContrast(MRIaggr.Pat1_red,param="MASK_DWI_t0")
coords <- selectCoords(MRIaggr.Pat1_red)

#### 1- compute spatial groups using coordinates ####
res3DN18 <- calcGroupsCoords(coords=coords[MASK_DWI_t0==1,],Neighborhood="3D_N18")
res3DN18$group_size

## display the lesion within the brain
multiplot(coords,contrast=MASK_DWI_t0,legend=FALSE,
          index1=list(coords=coords,outline=TRUE,cex=0.5),
           num.main=FALSE)

## display the lesion spatial groups
multiplot(res3DN18$df.group[,c("i","j","k")],contrast=res3DN18$df.group[,"group"],
             legend=FALSE,cex=0.5,index1=list(coords=coords,outline=TRUE))

## zoom
multiplot(res3DN18$df.group[,c("i","j","k")],contrast=res3DN18$df.group[,"group"],
             legend=FALSE,num=1,cex=0.1,index1=list(coords=coords,outline=TRUE))

#### 2-compute spatial groups using an array ####
A.MASK_DWI_t0 <- df2array(MASK_DWI_t0,coords=coords)$contrast[[1]]
A.MASK_DWI_t0[A.MASK_DWI_t0==FALSE] <- NA

## display
graphics::image(A.MASK_DWI_t0[,,3])

## computation of the spatial groups
res3DN18.bis <- calcGroupsCoords(array=A.MASK_DWI_t0,Neighborhood="3D_N18")

res3DN18$group_size-res3DN18.bis$group_size # same result
}
\concept{calc.}
\keyword{functions}


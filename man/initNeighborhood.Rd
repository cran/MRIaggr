\name{initNeighborhood}
\title{Initialization of a neighborhood filter}
\alias{initNeighborhood}
\usage{
initNeighborhood(Neighborhood,method)
}
\description{
Return the neighborhood configuration corresponding to the specified name.
}
\arguments{
  \item{Neighborhood}{the name of neighborhood configuration. \cr Any of \code{"2D_N4"}, \code{"2D_N8"} \code{"3D_N6"} \code{"3D_N10"} \code{"3D_N18"} \code{"3D_N26"}.}
  \item{method}{the name of the function that called the initializer. \emph{character}.}
}
\details{
ARGUMENTS : \cr
 \code{Neighborhood} refers to classical neighborhood configurations : \cr
 The first two characters refer to the dimension d of the filter : \code{"2D"} or \code{"3D"}. \cr
 The third character must be \code{"_"}. \cr
 The fourth character refers to the type of filter and must be \code{"N"}. \cr
 The last one or two characters indicates the number of neighbors (denoted \code{n}) in each neighborhood.
}
\value{
 A n*d \emph{matrix} with in line the coordinates of the neighbors relative to the current observation.
}
\examples{
# 2D neighborhood
initNeighborhood("2D_N4",method="calcFilter") # rock neighborhood
initNeighborhood("2D_N8",method="calcFilter") # queen neighborhood

# 3D neighborhood
initNeighborhood("3D_N6",method="calcFilter") # rock neighborhood
initNeighborhood("3D_N10",method="calcFilter") # queen neighborhood
initNeighborhood("3D_N18",method="calcFilter") # queen neighborhood
initNeighborhood("3D_N26",method="calcFilter") # queen neighborhood

}
\concept{init.}
\keyword{functions}

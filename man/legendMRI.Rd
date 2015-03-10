\name{legendMRI}
\title{Diplay a legend of a contrast map}
\alias{legendMRI}
\description{
Display a legend from break values. For internal use.
}
\usage{
legendMRI(breaks,palette,mar,cex,cex.main,main,quantiles,digit)
}
\arguments{
  \item{breaks}{the break points to use to generate the color intervals. \emph{numeric vector}.}
  \item{palette}{the colors or the palette to use when associating colors to observation intensities. \emph{character vector}.}
 
  \item{mar}{the number of margin lines to be specified on the four sides of the legend. \emph{positive numeric vector of size 4}.}
  \item{cex.main}{the expansion factor for the legend title. \emph{positive numeric}.}
  \item{cex}{the expansion factor for the legend labels. \emph{positive numeric}.}
  \item{main}{an overall title for the legend. \emph{character}.}
  \item{quantiles}{the quantiles values to display on the legend. \emph{numeric vector of size 5} or \code{NULL}.}
  \item{digit}{the number of decimal places to use for the legend label. \emph{positive integer}.}
}
\details{
ARGUMENTS : \cr
Information about the \code{mar} argument can be found in \code{\link{par}}.

Information about the \code{palette} argument can be found in the details section of \code{\link{initCol}}.
}
\concept{init.}
\keyword{function,internal}


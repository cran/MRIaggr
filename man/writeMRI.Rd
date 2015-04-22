\name{writeMRI}
\title{Write an image file}
\alias{writeMRI}
\description{
Write an image file from an array.
}
\usage{
writeMRI(data, file, format, gzipped=TRUE, verbose=FALSE, size="NA_integer_")
}
\arguments{
  \item{data}{the data to export. \emph{three dimensional array}. REQUIRED.}
  \item{file}{the name of the image file. \emph{character}. REQUIRED.}
  \item{format}{the format of the image file. Can be \code{"raw.gz"}, \code{"analyze"},  \code{"nifti"} or \code{"dicom"}. REQUIRED.} 
  \item{size}{the number of bytes per element in the byte stream. \emph{integer}. See the documentation of the \code{writeBin} function for more details.}
  \item{gzipped}{is a character string that enables exportation of compressed (.gz) files.}
  \item{verbose}{is a logical variable that allows text-based feedback during execution of the function.} 
  }
\details{

This function requires to have installed the \emph{oro.nifti} package to work if argument format is set to analyze or nifti.

ARGUMENTS :
\code{file} argument corresponds to:
\itemize{
\item the \code{con} argument of the \code{base::readBin} function. 
\item the \code{filename} argument of the \code{oro.nifti::readANALYZE} function. It should be a pathname to .img or .hdr files without the suffix.
\item the \code{filename} argument of the \code{oro.nifti::readNIfTI} function. 
}

FUNCTION : \cr
This function is inspired from the \code{writeMRI} function of the mritc package. It calls
\code{writeANALYZE} of the oro.nifti package to write analyze files, 
\code{writeNIfTI} of the oro.nifti package to write nifti files, 
and \code{writeBin} of the base package to write rawb.gz.
}
\value{
  None.
}
\examples{ 
set.seed(10)
X <- rnorm(1000)

\dontrun{
writeMRI(array(X,c(10,10,10)),file="test",format="nifti",verbose=TRUE)
# range(array(X,c(10,10,10))-readMRI(file="test",format="nifti")[])

writeMRI(array(X,c(10,10,10)),file="test",format="analyze",verbose=TRUE)
# range(array(X,c(10,10,10))-readMRI(file="test",format="analyze")[])

writeMRI(array(X,c(10,10,10)),file="test",format="rawb.gz",verbose=TRUE,size="NA_integer_")
# range(array(X,c(10,10,10))-readMRI(file="test",
#       format="rawb.gz",size="NA_integer_",dim=c(10,10,10)))
}

}
\concept{convert.}
\keyword{functions}
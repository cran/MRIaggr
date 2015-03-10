\name{readMRI}
\title{Read imaging file}
\alias{readMRI}
\description{
Read an imaging file and convert it into an array.
}
\usage{
readMRI(file, format, what="integer", na.value = 0, dim = NULL,  
        SPM=FALSE, reorient = FALSE, flipud=FALSE)
}
\arguments{
  \item{file}{the file name of the imaging file. \emph{character}. REQUIRED.}
  \item{format}{the format of the image file. Can be \code{"raw.gz"}, \code{"analyze"},  \code{"nifti"} or \code{"dicom"}. REQUIRED.} 
  \item{what}{an object whose mode will give the mode of the vector to be read, or a character vector of length one describing the mode: one of "numeric", "double", "integer", "int", "logical", "complex", "character", "raw". Only active if \code{format} equals \code{rawb.gz}.}
  \item{na.value}{the value with which NA values are replaced. \emph{numeric} or NA.}    
  \item{dim}{the number of bytes per element in the byte stream. The default, NA_integer_, uses the natural size. Size changing is not supported for raw and complex vectors. Only active if \code{format} equals \code{rawb.gz}.}
  \item{SPM}{is a logical variable (default = FALSE) that forces the voxel data values to be rescaled using the funused1 ANALYZE header field. This is an undocumented convention of ANALYZE files processed using the Statistical Parametric Mapping (SPM) software. Only active if \code{format} equals \code{analyse}.}
  \item{reorient}{is a logical variable (default = TRUE) that enforces Qform/Sform transformations. Only active if \code{format} equals \code{nifti}.}
  \item{flipud}{is a logical variable for vertical flipping of the image (default is TRUE). Only active if \code{format} equals \code{dicom}.}
  }
\details{
ARGUMENTS :
\code{file} argument corresponds to:
\itemize{
\item the \code{con} argument of the \code{base::readBin} function. 
\item the \code{fname} argument of the \code{oro.nifti::readANALYZE} function. It should be a pathname to .img or .hdr files without the suffix.
\item the \code{fname} argument of the \code{oro.nifti::readNIfTI} function. 
\item the \code{fname} argument of the \code{oro.dicom::readDICOMFile} function. It should be the file name with suffix. 
}
The \code{what} and \code{dim} correspond to the \code{what} and \code{size}  argument of the \code{base::readBin} function. \cr
The \code{SPM} corresponds to the \code{SPM} argument of the \code{oro.nifti::readANALYZE} function. \cr
The \code{reorient} corresponds to the \code{reorient} argument of the \code{oro.nifti::readNIfTI} function. \cr
The \code{flipud} corresponds to the \code{flipud} argument of the \code{oro.dicom::readDICOMFile} function. \cr

FUNCTION : \cr
This function is a slightly modified version of the \code{readMRI} function from the mritc package.
}
\value{
  An \code{array}.
}
\examples{ 
## load a nifti file
path <- system.file(file.path("nifti"),package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path,"TTP_t0"),format="nifti")
dim(nifti.Pat1_TTP_t0)

\dontrun{
## load an analyse file (example of oro.nifti::readANALYZE)
path <- system.file("anlz", package="oro.nifti")
analyse.avg <- readMRI(file.path(path,"avg152T1"),format="analyze")
graphics::image(analyse.avg[,,45,1])

### load a nifti file (example of oro.nifti::readNIfTI)
nifti.ffd <- readMRI(file.path(system.file("nifti", package="oro.nifti"),"filtered_func_data"),
                     format="nifti")
graphics::image(nifti.ffd[,,10,32])

## load a dicom file (examples of oro.dicom::readDICOMFile)
dicom.Abdo <- readMRI(system.file("dcm/Abdo.dcm", package="oro.dicom"),format="dicom")
graphics::image(dicom.Abdo, col=grey(0:64/64), axes=FALSE, xlab="", ylab="",
                main="Abdo.dcm")
}
}
\concept{convert.}
\keyword{functions}
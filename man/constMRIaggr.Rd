\name{constMRIaggr}
\title{Array constructor for MRIaggr object}
\alias{constMRIaggr}
\description{
Construct a \code{\linkS4class{MRIaggr}} object from a list of array, each array corresponding to a different contrast parameters.
}
\usage{
constMRIaggr(ls.array,identifier,param,default_value=NULL,
            pos_default_value=c(1,1,1),tol=10^{-10},trace=TRUE,rm.ls.array=FALSE) 
}
\arguments{
  \item{ls.array}{the value of the contrast parameter(s) for each observation. \code{list of array}. REQUIRED.}
  \item{identifier}{the identifier of the patient to which belong the contrast parameters. \emph{character}. REQUIRED.}
  \item{param}{the contrast parameter(s). \emph{character vector} or \code{NULL}. REQUIRED.}
  \item{default_value}{the reference values of the contrast parameters (e.g. the background value). \emph{character} or \code{NULL} leading to search the reference value in \code{array[pos_default_value]}.}
  \item{pos_default_value}{the coordinates of the observations that contains the reference value. \emph{numeric vector}.}   
  
  \item{tol}{numeric precision for the consistency check. \emph{positive numeric}.}
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
  \item{rm.ls.array}{should the object on which the \code{ls.array} argument points be removed form the global environment ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
All the array in \code{ls.array} in must have the same dimensions.

Information about the \code{param} argument can be found in the details section of \code{\link{initParameter}}.

\code{pos_default_value} is active only if \code{default_value} is set to \code{NULL}.
}
\value{
  A \code{\linkS4class{MRIaggr}} object.
}
\examples{
#### 1- 1st method ####
## load nifti files
path <- system.file(file.path("nifti"),package = "MRIaggr")

nifti.Pat1_TTP_t0 <- readMRI(file=file.path(path,"TTP_t0"),format="nifti")
nifti.Pat1_DWI_t0 <- readMRI(file=file.path(path,"DWI_t0"),format="nifti")
nifti.Pat1_MASK_DWI_t0 <- readMRI(file=file.path(path,"MASK_DWI_t0"),format="nifti")
nifti.Pat1_MASK_T2_FLAIR_t2 <- readMRI(file=file.path(path,"MASK_T2_FLAIR_t2"),format="nifti")

## convert them to MRIaggr
MRIaggr.Pat1 <- constMRIaggr(list(nifti.Pat1_TTP_t0,nifti.Pat1_DWI_t0,
                                  nifti.Pat1_MASK_DWI_t0,nifti.Pat1_MASK_T2_FLAIR_t2),
                identifier="Pat1",param=c("TTP_t0","DWI_t0","MASK_DWI_t0","MASK_T2_FLAIR_t2"))


#### 2- 2nd method ####
## load nifti files
param <- c("DWI_t0.nii","MASK_DWI_t0.nii","MTT_t0.nii","TTP_t0.nii","T1_t0.nii","T2_GRE_t0.nii",
           "MTT_t1.nii","TTP_t1.nii","T2_FLAIR_t2.nii","MASK_T2_FLAIR_t2.nii")

ls.array <- list()
for(iter_param in 1:length(param)){
  ls.array[[iter_param]] <- readMRI(file=file.path(path,param[iter_param]),format="nifti")
}

## convert them to MRIaggr
param <- gsub(".nii","",param)

MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1",param=param)

#### additionnal examples
\dontrun{
## load an analyse file (example of oro.nifti::readANALYZE)
path <- system.file("anlz", package="oro.nifti")
analyse.avg <- readMRI(file.path(path,"avg152T1"),format="analyze")
analyse.MRIaggr <- constMRIaggr(analyse.avg,param="avg",identifier="Pat1")

### load a nifti file (example of oro.nifti::readNIfTI)
nifti.ffd <- readMRI(file.path(system.file("nifti", package="oro.nifti"),"filtered_func_data"),
                     format="nifti")
nifti.MRIaggr <- constMRIaggr(nifti.ffd[,,,32,drop=FALSE],param="ffd",identifier="Pat1")

## load a dicom file (examples of oro.dicom::readDICOMFile)
dicom.Abdo <- readMRI(system.file("dcm/Abdo.dcm", package="oro.dicom"),format="dicom")
dicom.MRIaggr <- constMRIaggr(dicom.Abdo,param="Abdo",identifier="Pat1")
}
}
\concept{const.}
\keyword{functions}


\docType{class}
\name{MRIaggr-class}
\title{Class "MRIaggr"}
\alias{MRIaggr}
\alias{MRIaggr-class}
\description{
  Patient-specific storage of multivariate contrast data and clinical data.
}
\arguments{
  \item{identifier}{the patient identifier. \emph{character}.}  
  \item{contrast}{the value of the contrast parameters for each observation. \emph{data.frame}.}  
  \item{clinic}{the clinical data of the patient. \emph{data.frame}.}  
  \item{voxelDim}{the spatial dimensions of the lattice containing the observations. \emph{data.frame}.}  
  \item{default_value}{the reference values of the contrast parameters. \emph{data.frame}.}  

  \item{history}{the list of the \code{calc} or \code{const} methods that have been already applied on the \code{MRIaggr} object. \emph{data.frame}.}  
  \item{normalization}{the normalization values for the contrast parameters. \emph{list}.}  
  \item{hemispheres}{the presence or absence of lesion in each cerebral hemisphere. \emph{data.frame}.}  
  \item{midplane}{the position of the mid-sagittal plan. \emph{data.frame}.}  
  
  \item{table_lesion}{the vertical distribution of the lesion volumes. \emph{data.frame}.}  
  \item{table_reperfusion}{the volumic reperfusion data. \emph{data.frame}.}  
  \item{table_hypoperfusion}{the volumic hypoperfusion data. \emph{data.frame}.}  

  \item{ls_statDesc}{a slot to store additional data. \emph{list}.}  
}
\section{S4 methods}{
Getters : 
\describe{
 \item{\code{\link{selectContrast}}}{}
 \item{\code{\link{selectCoords}}}{}
 \item{\code{\link{selectClinic}}}{}
 \item{\code{\link{selectDefault_value}}}{}
 \item{\code{\link{selectDescStats}}}{}
 \item{\code{\link{selectHemispheres}}}{}
 \item{\code{\link{selectHistory}}}{}
 \item{\code{\link{selectIdentifier}}}{}
 \item{\code{\link{selectMidplane}}}{}
 \item{\code{\link{selectN}}}{}
 \item{\code{\link{selectNormalization}}}{}
 \item{\code{\link{selectParameter}}}{}
 \item{\code{\link{selectTable}}}{}
 \item{\code{\link{selectVoxelDim}}}{}
 \item{\code{\link{selectVoxelSize}}}{}
}

Setters : 
\describe{
\item{\code{\link{affectContrast<-}}}{}
\item{\code{\link{affectClinic<-}}}{}
\item{\code{\link{affectDescStats<-}}}{}
\item{\code{\link{affectHemisphere<-}}}{}
\item{\code{\link{affectNormalization<-}}}{}
\item{\code{\link{affectTable<-}}}{}
\item{\code{\link{supprContrast<-}}}{}
\item{\code{\link{supprDescStats<-}}}{}
}

Calculators : 
\describe{
\item{\code{\link{calcBrainMask}}}{}
\item{\code{\link{calcControlateral}}}{}
\item{\code{\link{calcDistMask}}}{}
\item{\code{\link{calcDistTissues}}}{}
\item{\code{\link{calcFilter}}}{}
\item{\code{\link{calcGroupsMask}}}{}
\item{\code{\link{calcHemisphere}}}{}
\item{\code{\link{calcROCthreshold}}}{}
\item{\code{\link{calcNormalization}}}{}
\item{\code{\link{calcRegionalContrast}}}{}
\item{\code{\link{calcSmoothMask}}}{}
\item{\code{\link{calcTableHypoReperf}}}{}
\item{\code{\link{calcTableLesion}}}{}
\item{\code{\link{calcThresholdMRIaggr}}}{}
\item{\code{\link{calcTissueType}}}{}
\item{\code{\link{calcW}}}{}
}

Displayers : 
\describe{
\item{\code{\link{boxplotMask}}}{}
\item{\code{\link{heatmapMRIaggr}}}{}
\item{\code{\link{multiplot}}}{}
\item{\code{\link{plotDistClass}}}{}
\item{\code{\link{pointsHemisphere}}}{}
\item{\code{\link{plotLesion3D}}}{}
\item{\code{\link{plotTableLesion}}}{}
\item{\code{\link{outlineMRIaggr}}}{}
\item{\code{\link{summary,MRIaggr-method}}}{}
}

Constructors : 
\describe{
\item{\code{\link{constCompressMRIaggr}}}{}
\item{\code{\link{constReduceMRIaggr}}}{}
}

Initializers : 
\describe{
\item{\code{\link{initNum}}}{} 
\item{\code{\link{initParameter}}}{}  
}
}
\seealso{
\code{\link{constMRIaggr}} to build a \code{\linkS4class{MRIaggr}} object from a list of array.  \cr
\code{\link{readMRI}} to read imaging files. 
}
\keyword{classes}

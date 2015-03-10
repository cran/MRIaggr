\name{calcContro_cpp}
\title{Compute controlateral normalization}
\alias{calcContro_cpp}
\usage{
calcContro_cpp(contrast, coords_px, index_k, index_k_contro, d_lim, 
         lambda, param_ref, var_ref, type_moy, type_med, type_NN, trace)
}
\description{
C++ function called by \code{\link{calcControlateral}} to seek the controlateral voxel and the controlateral values. For internal use.
}
\arguments{
  \item{contrast}{the contrast value of each voxel on a given slice. \emph{matrix}.}
  \item{coords_px}{the coordinates of the observations. \emph{matrix}.}
  \item{index_k}{the index of the observations on the hemisphere of interest. \emph{integer vector}.}
  \item{index_k_contro}{the index of the observations on the controlateral hemisphere. \emph{integer vector}.}
  \item{d_lim}{the distance within which the controlateral values are considered. \emph{numeric}.}
  \item{lambda}{the importance of the penalization by the distance. \emph{numeric}.}
  \item{param_ref}{the parameter to be used as a reference to identify the controlateral voxel. \emph{character}.}
  \item{var_ref}{the variance of the reference parameter. \emph{numeric}.}
  \item{type_moy}{should the mean controlateral value be used ? \emph{logical}.}
  \item{type_med}{should the median controlateral value be used ? \emph{logical}.}
  \item{type_NN}{should the closest controlateral voxel according to the reference parameter be used ? \emph{logical}.}
  \item{trace}{should the execution of the function be traced ? \emph{logical}.}
}
\concept{cpp.}
\keyword{functions,internal}

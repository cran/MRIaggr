\name{summary-methods}
\title{Summary Method for Class "MRIaggr"}
\alias{summary,MRIaggr-method}
\description{
 Summarize some information of an object of class "MRIaggr". 
}
\usage{
\S4method{summary}{MRIaggr}(object,param=FALSE,clinic=FALSE,descStats=FALSE,history=FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{should detailed information be printed for the contrast parameters ? \emph{logical}.}
  \item{clinic}{should detailed information be printed for the clinical attribute ? \emph{logical}.}
  \item{descStats}{should detailed information be printed for the \code{ls_descStats} attribute ? \emph{logical}.}
  \item{history}{should the\code{calc} and \code{const} methods that have been applied to the object be listed ? \emph{logical}.}
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

summary(MRIaggr.Pat1_red)
}
\concept{plot.}
\keyword{methods,summary}


\name{specify_sample}
\alias{specify_sample}
\title{
  create an object of class 'popObj' required for further processing
}
\description{
  adjust sampling weights to given totals based on household-level and/or individual level constraints
}
\usage{
specify_pop(data, hhid, strata=NULL)
}
\arguments{
  \item{data}{a \code{data.frame} or \code{data.table} featuring sample data.}
  \item{hhid}{character vector of length 1 specifying variable containing household ids within \code{samp}.}
  \item{strata}{character vector of length 1 specifing variable name within \code{samp} of variable holding information on strata, e.g. regions or NULL if such variable does not exist.} 
}

\examples{
data(eusilcP)
p <- specify_pop(data=eusilcP, hhid="db030", strata="db040")
print(class(p))
}
\keyword{method}

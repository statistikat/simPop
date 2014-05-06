\name{specify_sample}
\alias{specify_sample}
\title{
  create an object of class 'sampleObj' required for further processing
}
\description{
  adjust sampling weights to given totals based on household-level and/or individual level constraints
}
\usage{
specify_sample(data, hhid, weight, strata=NULL)
}
\arguments{
  \item{data}{a \code{data.frame} or \code{data.table} featuring sample data.}
  \item{hhid}{character vector of length 1 specifying variable containing household ids within \code{samp}.}
  \item{weight}{character vector of length 1 specifying variable holding sampling weights within \code{samp}.}
  \item{strata}{character vector of length 1 specifing variable name within \code{samp} of variable holding information on strata, e.g. regions or NULL if such variable does not exist.} 
}

\examples{
data(eusilcS)
s <- specify_sample(data=eusilcS, hhid="db030", weight="rb050", strata="db040")
print(class(s))
}
\keyword{method}

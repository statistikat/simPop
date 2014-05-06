\name{specify_sample}
\alias{specify_sample}
\title{
  create an object of class 'popObj' required for further processing
}
\description{
  create an standardized population object of class 'popObj' containing information on household ids, household sizes, person ids and optionally strata.
}
\usage{
specify_pop(data, hhid, hhsize=NULL, pid=NULL, strata=NULL)
}
\arguments{
  \item{data}{a \code{data.frame} or \code{data.table} featuring sample data.}
  \item{hhid}{character vector of length 1 specifying variable containing household ids within \code{samp}.}
  \item{hhsize}{character vector of length 1 specifying variable containing household sizes within \code{samp}. If NULL, household sizes are automatically calculated.}
  \item{pid}{character vector of length 1 specifying variable containing person ids within \code{samp}. If NULL, person ids are automatically calculated.}
  \item{strata}{character vector of length 1 specifing variable name within \code{samp} of variable holding information on strata, e.g. regions or NULL if such variable does not exist.} 
}

\examples{
data(eusilcP)
p <- specify_pop(data=eusilcP, hhid="db030", strata="db040")
print(class(p))
}
\keyword{method}

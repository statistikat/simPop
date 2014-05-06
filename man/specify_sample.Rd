\name{specify_sample}
\alias{specify_sample}
\title{
  create an object of class 'sampleObj' required for further processing
}
\author{Bernhard Meindl}
\description{
  create an standardized sample object of class 'sampObj' containing information on weights, household ids, household sizes, person ids and optionally strata.
}
\usage{
specify_sample(data, hhid, hhsize=NULL, pid=NULL, weight, strata=NULL, additional=NULL)
}
\arguments{
  \item{data}{a \code{data.frame} or \code{data.table} featuring sample data.}
  \item{hhid}{character vector of length 1 specifying variable containing household ids within \code{data}.}
  \item{hhsize}{character vector of length 1 specifying variable containing household sizes within \code{samp}. If NULL, household sizes are automatically calculated.}
  \item{pid}{character vector of length 1 specifying variable containing person ids within \code{data}. If NULL, person ids are automatically calculated.}
  \item{weight}{character vector of length 1 specifying variable holding sampling weights within \code{data}.}
  \item{strata}{character vector of length 1 specifing variable name within \code{data} of variable holding information on strata, e.g. regions or NULL if such variable does not exist.}
  \item{additional}{character vector of additional variables withing \code{data} that are important for the household structure, typically avariables such as 'age' or 'sex'.}
}
\examples{
data(eusilcS)
s <- specify_sample(data=eusilcS, hhid="db030", weight="rb050", strata="db040", additional=c("rb090", "age"))
print(class(s))
}
\keyword{method}

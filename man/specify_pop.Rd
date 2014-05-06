\name{specify_pop}
\alias{specify_pop}
\title{
  create an object of class 'popObj' required for further processing
}
\author{Bernhard Meindl}
\description{
  create an standardized population object of class 'popObj' containing information on household ids, household sizes, person ids and optionally strata.
}
\usage{
specify_pop(data, hhid, hhsize=NULL, pid=NULL, strata=NULL, additional=NULL)
}
\arguments{
  \item{data}{a \code{data.frame} or \code{data.table} featuring population data.}
  \item{hhid}{character vector of length 1 specifying variable containing household ids within \code{data}.}
  \item{hhsize}{character vector of length 1 specifying variable containing household sizes within \code{data}. If NULL, household sizes are automatically calculated.}
  \item{pid}{character vector of length 1 specifying variable containing person ids within \code{data}. If NULL, person ids are automatically calculated.}
  \item{strata}{character vector of length 1 specifing variable name within \code{data} of variable holding information on strata, e.g. regions or NULL if such variable does not exist.} 
  \item{additional}{character vector of additional variables withing \code{data} that are important for the household structure, typically avariables such as 'age' or 'sex'.}
}
\examples{
data(eusilcP)
p <- specify_pop(data=eusilcP, hhid="hid", strata="region", additional=c("age", "gender", "citizenship"))
print(class(p))
}
\keyword{method}

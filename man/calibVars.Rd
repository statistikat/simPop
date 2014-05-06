\name{calibVars}
\alias{calibVars}
\title{Construct a matrix of binary variables for calibration}
\usage{
  calibVars(x)
}
\arguments{
  \item{x}{a vector that can be interpreted as factor, or a
  matrix or \code{data.frame} consisting of such
  variables.}
}
\value{
  A matrix of binary variables that indicate membership to
  the corresponding factor levels.
}
\description{
  Construct a matrix of binary variables for calibration of
  sample weights according to known marginal population
  totals.
}
\examples{
data(eusilcS)
# default method
aux <- calibVars(eusilcS$rb090)
head(aux)
# data.frame method
aux <- calibVars(eusilcS[, c("db040", "rb090")])
head(aux)
}
\author{
  Bernhard Meindl and Andreas Alfons
}
\seealso{
  \code{\link{calibSample}}
}
\keyword{survey}


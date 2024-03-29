% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simStructure.R
\name{simStructure}
\alias{simStructure}
\title{Simulate the household structure of population data}
\usage{
simStructure(
  dataS,
  method = c("direct", "multinom", "distribution"),
  basicHHvars,
  seed = 1,
  MaxNWts = 1e+07
)
}
\arguments{
\item{dataS}{an object of class \code{dataObj} containing household survey
data that is usually generated with \code{\link{specifyInput}}.}

\item{method}{a character string specifying the method to be used for
simulating the household sizes.  Accepted values are \code{"direct"}
(estimation of the population totals for each combination of stratum and
household size using the Horvitz-Thompson estimator), \code{"multinom"}
(estimation of the conditional probabilities within the strata using a
multinomial log-linear model and random draws from the resulting
distributions), or \code{"distribution"} (random draws from the observed
conditional distributions within the strata).}

\item{basicHHvars}{a character vector specifying important variables for the
household structure that need to be available in \code{dataS}. Typically
variables such as age or sex may be used.}

\item{seed}{optional; an integer value to be used as the seed of the random
number generator, or an integer vector containing the state of the random
number generator to be restored.}

\item{MaxNWts}{optional; an integer value for the multinom method for controlling 
the maximum number of weights.}
}
\value{
An object of class \code{simPopObj} containing the simulated
population household structure as well as the underlying sample that was
provided as input.
}
\description{
Simulate basic categorical variables that define the household structure
(typically variables such as household ID, age and gender) of population
data by resampling from survey data.
}
\note{
The function \code{\link{sample}} is used, which gives results
incompatible with those from < 2.2.0 and produces a warning the first time
this happens in a session.
}
\examples{

data(eusilcS)
\dontrun{
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
eusilcP <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
class(eusilcP)
eusilcP
}

}
\references{
M. Templ, B. Meindl, A. Kowarik, A. Alfons, O. Dupriez (2017) Simulation of Synthetic Populations for Survey Data Considering Auxiliary
Information. \emph{Journal of Statistical Survey}, \strong{79} (10), 1--38. \doi{10.18637/jss.v079.i10}
}
\seealso{
\code{\link{simCategorical}}, \code{\link{simContinuous}},
\code{\link{simComponents}}, \code{\link{simEUSILC}}
}
\author{
Bernhard Meindl and Andreas Alfons
}
\keyword{datagen}

\name{spCdfplot}
\Rdversion{1.1}
\alias{spCdfplot}
\alias{spCdfplot.default}
\title{
  Plot weighted cumulative distribution functions
}
\description{
  Plot cumulative distribution functions, possibly broken down according to conditioning variables and taking into account sample weights.
}
\usage{
spCdfplot(inp, x, cond = NULL, approx = c(FALSE, TRUE), 
  n = 10000, bounds = TRUE, \dots)
}
\arguments{
  \item{inp}{an object of class \code{\linkS4class{synthPopObj}} containing survey sample and synthetic population data.}
  \item{x}{a character vector specifying the columns of data available in the sample and the population (specified in input object 'inp') to be plotted.}
  \item{cond}{an optional character vector (of length 1, if used) specifying the conditioning variable.}
  \item{approx}{logicals indicating whether approximations of the cumulative distribution functions should be computed.  The default is to use \code{FALSE} for the survey data and \code{TRUE} for the population data.}
  \item{n}{integers specifying the number of points at which the approximations take place (see \code{\link[stats:approxfun]{approx}}).  It is used wherever \code{approx} is \code{TRUE}.}
  \item{bounds}{a logical indicating whether vertical lines should be drawn at 0 and 1 (the bounds for cumulative distribution functions).}
  \item{\dots}{further arguments to be passed to \code{\link[lattice]{xyplot}}.}
}
\details{
  Weights are directly extracted from the input object \code{inp} and are taken into account by adjusting the step height.  To be precise, the weighted step height for an observation is defined as its weight divided by the sum of all weights\eqn{\ ( w_{i} / \sum_{j = 1}^{n} w_{j} ).}{.}
}
\value{
  An object of class \code{"trellis"}, as returned by \code{\link[lattice]{xyplot}}.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{spCdf}}, \code{\link[lattice]{xyplot}}
}
\examples{
\dontrun{

## these take some time and are not run automatically
## copy & paste to the R command line

set.seed(1234)  # for reproducibility
data(eusilcS)   # load sample data
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", 
  strata="db040", weight="db090")
synthPop <- simStructure(data=inp, method="direct",
  basicHHvars=c("age", "rb090", "hsize", "pl030", "pb220a"))

# multinomial model with random draws
eusilcM <- simContinuous(synthPop, additional="netIncome", 
  upper=200000, equidist=FALSE)
class(eusilcM)

# plot results
spCdfplot(eusilcM, "netIncome", cond=NULL)
spCdfplot(eusilcM, "netIncome", cond="rb090", layout=c(1,2))
}
}
\keyword{hplot}

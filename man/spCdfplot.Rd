\name{spCdfplot}
\Rdversion{1.1}
\alias{spCdfplot}
\alias{spCdfplot.default}
\title{
  Plot (weighted empirical) cumulative distribution functions
}
\description{
  Plot (weighted empirical) cumulative distribution functions for survey 
  and population data, possibly broken down according to conditioning 
  variables.  For survey data, sample weights are taken into account.
}
\usage{
spCdfplot(x, \dots)

\method{spCdfplot}{default}(x, weights = NULL, cond = NULL, dataS, dataP = NULL, 
          approx = NULL, n = 10000, bounds = TRUE, \dots)
}
\arguments{
  \item{x}{for the default method (currently the only method implemented), a 
    character vector specifying the columns of \code{dataS} and \code{dataP} to 
    be plotted.}
  \item{weights}{a character string specifying the column of \code{dataS} that 
    contains the (personal) sample weights.}
  \item{cond}{an optional character vector specifying conditioning variables.}
  \item{dataS}{a \code{data.frame} containing household survey data.}
  \item{dataP}{optional; a \code{data.frame} containing simulated population 
    data or a list of such \code{data.frame}s.}
  \item{approx}{logicals indicating whether approximations of the cumulative 
    distribution functions should be computed.  The default is to use 
    \code{FALSE} for the survey data and \code{TRUE} for any population data.  
    If no population data are supplied, a single logical should be used.  On 
    the other hand, if any population data are supplied, the behavior is as 
    follows.  If a single logical is supplied, it is used for the population 
    data and \code{FALSE} is used for the survey data.  If a vector of length 
    two is supplied, the first value is used for the survey data and the second 
    for the population data.  Note that if multiple populations are supplied, 
    the same value is used for all of them.}
  \item{n}{integers specifying the number of points at which the approximations 
    for the respective data sets take place (see 
    \code{\link[stats:approxfun]{approx}}).  If a single value is supplied, it 
    is used wherever \code{approx} is \code{TRUE}.  If a vector of length two 
    and any population data are supplied, the first value is used for the 
    survey data and the second for the population data (in case the 
    corresponding values of \code{approx} are \code{TRUE}).  Note that if 
    multiple populations are supplied, the same value is used for all of them.}
  \item{bounds}{a logical indicating whether vertical lines should be drawn at 
    0 and 1 (the bounds for cumulative distribution functions).}
  \item{\dots}{for the generic function, further arguments to be passed down to 
    methods.  For the default method, further arguments to be passed to 
    \code{\link[lattice]{xyplot}}.}
}
\details{
  Sample weights are taken into account by adjusting the step height.  To be 
  precise, the weighted step height for an observation is defined as its weight 
  divided by the sum of all weights\eqn{\ ( w_{i} / \sum_{j = 1}^{n} w_{j} ).}{.}
}
\value{
  An object of class \code{"trellis"}, as returned by 
  \code{\link[lattice]{xyplot}}.
}
\author{
  Andreas Alfons
}
\note{
  A formula interface may be added in the future.
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
# multinomial model with random draws
eusilcM <- simEUSILC(eusilcS, upper = 200000, equidist = FALSE)
# two-step regression
eusilcT <- simEUSILC(eusilcS, method = "twostep")
# plot results by gender
spCdfplot("netIncome", "rb050", "rb090", dataS = eusilcS, 
    dataP = list(M = eusilcM, T = eusilcT), layout = c(1, 2))
}
}
\keyword{hplot}

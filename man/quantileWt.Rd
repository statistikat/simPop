\name{quantileWt}
\alias{quantileWt}
\alias{quantileWt.default}
\alias{quantileWt.dataObj}
\title{
  Weighted sample quantiles
}
\description{
  Compute quantiles taking into account sample weights.
}
\usage{
quantileWt(x, \dots)

\method{quantileWt}{default}(x, weights=NULL, probs=seq(0, 1, 0.25), na.rm=TRUE, \dots)

\method{quantileWt}{dataObj}(x, vars, probs=seq(0, 1, 0.25), na.rm = TRUE, \dots)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{weights}{an optional numeric vector containing sample weights.}
  \item{vars}{a character vector of length 1 specifying a variable name that is available in the data-slot of \code{x} and which is used for the calculation.}
  \item{probs}{a numeric vector of probabilities with values in \eqn{[0, 1]}.}
  \item{na.rm}{a logical indicating whether any \code{NA} or \code{NaN} values
    should be removed from \code{x} before the quantiles are computed.  Note
    that the default is \code{TRUE}, contrary to the function
    \code{\link[stats]{quantile}}.}
  \item{\dots}{for the generic function \code{quantileWt} additional arguments 
    to be passed to methods.  Additional arguments not included in the definition 
    of the methods are currently ignored.}    
}
\details{
  If weights are not specified then \code{quantile(x, probs, na.rm=na.rm,
  names=FALSE, type=1)} is used for the computation.

  Note probabilities outside \eqn{[0, 1]} cause an error.
}
\value{
  A vector of the (weighted) sample quantiles.
}
\author{
  Stefan Kraft and Bernhard Meindl

  A basic version of this function was provided by Cedric Beguin and Beat
  Hulliger.
}
\seealso{\code{\link[stats]{quantile}}}
\examples{
data(eusilcS)
(quantileWt(eusilcS$netIncome, weights=eusilcS$rb050))

# dataObj-method
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
(quantileWt(inp, vars="netIncome"))
}
\keyword{univar}

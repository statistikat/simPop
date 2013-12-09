\name{quantileWt}
\Rdversion{1.1}
\alias{quantileWt}
\title{
  Weighted sample quantiles
}
\description{
  Compute quantiles taking into account sample weights.
}
\usage{
quantileWt(x, weights = NULL, probs = seq(0, 1, 0.25), na.rm = TRUE)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{weights}{an optional numeric vector containing sample weights.}
  \item{probs}{a numeric vector of probabilities with values in \eqn{[0, 1]}.}
  \item{na.rm}{a logical indicating whether any \code{NA} or \code{NaN} values 
    should be removed from \code{x} before the quantiles are computed.  Note 
    that the default is \code{TRUE}, contrary to the function 
    \code{\link[stats]{quantile}}.}
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
  Stefan Kraft

  A basic version of this function was provided by Cedric Beguin and Beat 
  Hulliger.
}
\seealso{\code{\link[stats]{quantile}}}
\examples{
data(eusilcS)
quantileWt(eusilcS$netIncome, weights=eusilcS$rb050)
}
\keyword{univar}

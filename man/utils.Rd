\name{utils}
\Rdversion{1.1}
\alias{meanWt}
\alias{varWt}
\alias{covWt}
\alias{covWt.default}
\alias{covWt.matrix}
\alias{covWt.data.frame}
\alias{corWt}
\alias{corWt.default}
\alias{corWt.matrix}
\alias{corWt.data.frame}
\title{
  Weighted mean, variance, covariance matrix and correlation matrix
}
\description{
  Compute mean, variance, covariance matrix and correlation matrix, taking into 
  account sample weights.
}
\usage{
meanWt(x, weights, na.rm = TRUE)

varWt(x, weights, na.rm = TRUE)

covWt(x, \dots)

\method{covWt}{default}(x, y, weights, \dots)

\method{covWt}{matrix}(x, weights, \dots)

\method{covWt}{data.frame}(x, weights, \dots)

corWt(x, \dots)

\method{corWt}{default}(x, y, weights, \dots)

\method{corWt}{matrix}(x, weights, \dots)

\method{corWt}{data.frame}(x, weights, \dots)
}
\arguments{
  \item{x}{for \code{meanWt} and \code{varWt}, a numeric vector.  For 
    \code{covWt} and \code{corWt}, a numeric vector, matrix or 
    \code{data.frame}.}
  \item{y}{a numeric vector.  If missing, this defaults to \code{x}.}
  \item{weights}{an optional numeric vector containing sample weights.}
  \item{na.rm}{a logical indicating whether any \code{NA} or \code{NaN} values 
    should be removed from \code{x} before computation.  Note that the default 
    is \code{TRUE}.}
  \item{\dots}{for the generic functions \code{covWt} and \code{corWt}, 
    additional arguments to be passed to methods.  Additional arguments not 
    included in the definition of the methods are ignored.}
}
\details{
  \code{meanWt} is a simple wrapper that calls \code{mean(x, na.rm=na.rm)} if 
  \code{weights} is missing and \code{weighted.mean(x, w=weights, na.rm=na.rm)} 
  otherwise.
  
  \code{varWt} calls \code{var(x, na.rm=na.rm)} if \code{weights} is missing.
  
  \code{covWt} and \code{corWt} always remove missing values pairwise and 
  call \code{cov} and \code{cor}, respectively, if \code{weights} is missing.
}
\value{
  For \code{meanWt}, the (weighted) mean.
  
  For \code{varWt}, the (weighted) variance.
  
  For \code{covWt}, the (weighted) covariance matrix or, for the default 
  method, the (weighted) covariance.
  
  For \code{corWt}, the (weighted) correlation matrix or, for the default 
  method, the (weighted) correlation coefficient.
}
\author{
  Stefan Kraft and Andreas Alfons
}
\seealso{
  \code{\link{mean}}, \code{\link[stats]{weighted.mean}}, 
  \code{\link[stats:cor]{var}}, \code{\link[stats:cor]{cov}}, 
  \code{\link[stats:cor]{cor}}
}
\examples{
data(eusilcS)
meanWt(eusilcS$netIncome, weights=eusilcS$rb050)
sqrt(varWt(eusilcS$netIncome, weights=eusilcS$rb050))
}
\keyword{univar}
\keyword{multivariate}
\keyword{array}
\encoding{utf8}
\name{calibWeights}
\alias{calibWeights}
\title{Calibrate sample weights}
\usage{
  calibWeights(X, d, totals, q = NULL,
    method = c("raking", "linear", "logit"),
    bounds = c(0, 10), maxit = 500, tol = 1e-06,
    eps = .Machine$double.eps)
}
\arguments{
  \item{X}{a matrix of binary calibration variables (see
  \code{\link{calibVars}}).}

  \item{d}{a numeric vector giving the initial sample
  weights.}

  \item{totals}{a numeric vector of population totals
  corresponding to the calibration variables in \code{X}.}

  \item{q}{a numeric vector of positive values accounting
  for heteroscedasticity.  Small values reduce the
  variation of the \emph{g}-weights.}

  \item{method}{a character string specifying the
  calibration method to be used.  Possible values are
  \code{"linear"} for the linear method, \code{"raking"}
  for the multiplicative method known as raking and
  \code{"logit"} for the logit method.}

  \item{bounds}{a numeric vector of length two giving
  bounds for the g-weights to be used in the logit method.
  The first value gives the lower bound (which must be
  smaller than or equal to 1) and the second value gives
  the upper bound (which must be larger than or equal to
  1).}

  \item{maxit}{a numeric value giving the maximum number of
  iterations.}

  \item{tol}{the desired accuracy for the iterative
  procedure.}

  \item{eps}{the desired accuracy for computing the
  Moore-Penrose generalized inverse (see
  \code{\link[MASS]{ginv}}).}
}
\value{
  A numeric vector containing the \emph{g}-weights.
}
\description{
  Calibrate sample weights according to known marginal
  population totals. Based on initial sample weights, the
  so-called \emph{g}-weights are computed by generalized
  raking procedures.
}
\details{
  The final sample weights need to be computed by
  multiplying the resulting \emph{g}-weights with the
  initial sample weights.
}
\note{
  This is a faster implementation of parts of
  \code{\link[sampling]{calib}} from package
  \code{sampling}.  Note that the default calibration
  method is raking and that the truncated linear method is
  not yet implemented.
}
\examples{
data(eusilcS)
# construct auxiliary 0/1 variables for genders
aux <- calibVars(eusilcS$rb090)
# population totals
totals <- c(3990798, 4191431)
# compute g-weights
g <- calibWeights(aux, eusilcS$rb050, totals)
# compute final weights
weights <- g * eusilcS$rb050
summary(weights)
}
\author{
  Andreas Alfons
}
\references{
  Deville, J.-C. and \enc{Särndal}{Saerndal}, C.-E. (1992)
  Calibration estimators in survey sampling. \emph{Journal
  of the American Statistical Association}, \bold{87}(418),
  376--382.

  Deville, J.-C., \enc{Särndal}{Saerndal}, C.-E. and
  Sautory, O. (1993) Generalized raking procedures in
  survey sampling. \emph{Journal of the American
  Statistical Association}, \bold{88}(423), 1013--1020.
}
\seealso{
  \code{\link{calibVars}}
}
\keyword{survey}


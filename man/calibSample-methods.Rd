\name{calibSample-methods}
\docType{methods}
\alias{calibSample}
\alias{calibSample-methods}
\alias{calibSample,dataObj_or_synthPopObj,list-method}
\alias{calibSample,data.frame,list-method}
\title{Calibrate sample weights}
\description{
  Calibrate sample weights according to known marginal population totals. Based on initial sample weights, the so-called \emph{g}-weights are computed by generalized raking procedures.
}
\section{Methods}{
\describe{The function provides methods with the following signatures.
  \item{\code{signature(inp="dataObj_or_synthPopObj", totals="list",...)}}{
  Argument 'inp' must be an object of class \code{\linkS4class{dataObj}} or \code{\linkS4class{synthPopObj}}.
  }
  \item{\code{signature(inp="data.frame", totals="list",...)}}{
  Argument 'inp' must be an object of class 'data.frame'.
  }
}}

\details{
  The methods return a list containing both the \emph{g}-weights as well as the final weights (initial sampling weights adjusted by the \emph{g}-weights.
}
\note{
  This is a faster implementation of parts of \code{\link[sampling]{calib}} from package \code{sampling}. Note that the default calibration method is raking and that the truncated linear method is not yet implemented.
}
\author{
  Andreas Alfons and Bernhard Meindl
}
\references{
  Deville, J.-C. and \enc{Saerndal}{Saerndal}, C.-E. (1992)
  Calibration estimators in survey sampling. \emph{Journal
  of the American Statistical Association}, \bold{87}(418),
  376--382.
  Deville, J.-C., \enc{Saerndal}{Saerndal}, C.-E. and
  Sautory, O. (1993) Generalized raking procedures in
  survey sampling. \emph{Journal of the American
  Statistical Association}, \bold{88}(423), 1013--1020.
}
\keyword{survey}

\examples{
data(eusilcS)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
eusilcP <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))

# synthPopObj-method
totals <- list()
totals[[1]] <- c(3990798, 4191431)
names(totals) <- "rb090"
res1 <- calibSample(inp=eusilcP, totals)

# data.frame-method
df <- as.data.frame(eusilcP@sample@data)
res2 <- calibSample(inp=df, totals, w="db090")

identical(res1, res2)
}
\keyword{methods}


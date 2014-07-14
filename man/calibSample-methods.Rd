\name{calibSample-methods}
\docType{methods}
\alias{calibSample}
\alias{calibSample-methods}
\alias{calibSample,synthPopObj,list-method}
\alias{calibSample,data.frame,list-method}
\title{Calibrate sample weights}
\description{
  Calibrate sample weights according to known marginal population totals. Based on initial sample weights, the so-called \emph{g}-weights are computed by generalized raking procedures.
}
\section{Methods}{
\describe{blabla
  \item{\code{signature(inp="synthPopObj", totals="list",...)}}{
  bla1
  }
  \item{\code{signature(inp="data.frame", totals="list",...)}}{
  bla2
  }
}}

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

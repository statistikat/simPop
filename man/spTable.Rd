\name{spTable}
\alias{spTable}
\title{
  Cross tabulations of expected and realized population sizes.
}
\description{
  Compute contingency tables of expected (i.e., estimated) and realized (i.e., simulated) population sizes. The expected values are obtained with the Horvitz-Thompson estimator.
}
\usage{
spTable(inp, select)
}
\arguments{
  \item{inp}{an object of class \code{\linkS4class{synthPopObj}} containing household survey and simulated population data.}
  \item{select}{character; vector defining the columns in slots 'pop' and 'sample' of argument 'input' that should be used for tabulation.}
}
\details{
  The contingency tables are computed with \code{\link{tableWt}}.
}
\value{
  A list of class \code{"spTable"} with the following components:
  \item{expected}{the contingency table estimated from the survey data.}
  \item{realized}{the contingency table computed from the simulated population
    data.}
}
\author{
  Andreas Alfons and Bernhard Meindl
}
\note{
  Sampling weights are automatically used from the input object 'inp'!
}
\seealso{
  \code{\link{spMosaic}}, \code{\link{tableWt}}
}
\examples{
data(eusilcS)
samp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize",
  strata="db040", weight="db090")
eusilcP <- simStructure(data=samp, method="direct", basicHHvars=c("age", "rb090"))
res <- spTable(eusilcP, select = c("age", "rb090"))
class(res)
}
\keyword{dplot}

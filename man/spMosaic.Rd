\name{spMosaic}
\alias{spMosaic}
\title{
  Mosaic plots of expected and realized population sizes
}
\description{
  Create mosaic plots of expected (i.e., estimated) and realized (i.e., simulated) population sizes.
}
\usage{
spMosaic(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{"spTable"} created using function \code{\link{spTable}}.}
  \item{\dots}{further arguments to be passed to \code{\link[vcd]{cotabplot}}.}
}
\details{
  The two tables of expected and realized population sizes are combined into a single table, with an additional conditioning variable indicating expected and realized values. A conditional plot of this table is then produced using
  \code{\link[vcd]{cotabplot}}.
}
\author{
  Andreas Alfons and Bernhard Meindl
}

\seealso{
  \code{\link{spTable}}, \code{\link[vcd]{cotabplot}}
}
\examples{
data(eusilcS)
samp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize",
  strata="db040", weight="db090")
eusilcP <- simStructure(data=samp, method="direct", basicHHvars=c("age","rb090"))
abb <- c("B","LA","Vi","C","St","UA","Sa","T","Vo")
tab <- spTable(eusilcP, select=c("rb090", "db040", "hsize"))
spMosaic(tab, labeling=labeling_border(abbreviate=c(db040=TRUE)))
}
\keyword{hplot}

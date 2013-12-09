\name{spMosaic}
\Rdversion{1.1}
\alias{spMosaic}
\alias{spMosaic.default}
\title{
  Mosaic plots of expected and realized population sizes
}
\description{
  Create mosaic plots of expected (i.e., estimated) and realized (i.e., 
  simulated) population sizes.
}
\usage{
spMosaic(x, \dots)

\method{spMosaic}{default}(x, weights = NULL, dataS, dataP, \dots)
}
\arguments{
  \item{x}{for the default method, an optional character vector specifying the 
    columns of \code{dataS} and \code{dataP} to be cross tabulated and plotted.  
    Otherwise, an object of class \code{"spTable"}.}
  \item{weights}{either a numeric vector containing the (personal) sample 
    weights, or a character string specifying the corresponding column of 
    \code{dataS}.}
  \item{dataS}{a \code{data.frame} containing household survey data.}
  \item{dataP}{a \code{data.frame} containing simulated population data.}
  \item{\dots}{further arguments to be passed to \code{\link[vcd]{cotabplot}}.}
}
\details{
  The two tables of expected and realized population sizes are combined into a 
  single table, with an additional contitioning variable indicating expected 
  and realized values.  A conditional plot of this table is then produced using 
  \code{\link[vcd]{cotabplot}}.
}
\author{
  Andreas Alfons
}
\note{
  A formula interface may be added in the future.
}
\seealso{
  \code{\link{spTable}}, \code{\link[vcd]{cotabplot}}
}
\examples{
set.seed(1234)  # for reproducibility
data(eusilcS)   # load sample data
eusilcP <- simStructure(eusilcS)
abb <- c("B","LA","Vi","C","St","UA","Sa","T","Vo")
spMosaic(c("rb090", "db040", "hsize"), "rb050", eusilcS, eusilcP, 
    labeling=labeling_border(abbreviate=c(db040=TRUE)))
}
\keyword{hplot}

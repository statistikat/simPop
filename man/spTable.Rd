\name{spTable}
\Rdversion{1.1}
\alias{spTable}
\alias{spTable.default}
\alias{spTable.formula}
\title{
  Cross tabulations of expected and realized population sizes
}
\description{
  Compute contingency tables of expected (i.e., estimated) and realized (i.e., 
  simulated) population sizes.  The expected values are obtained with the 
  Horvitz-Thompson estimator.
}
\usage{
spTable(dataS, dataP, select = NULL, weights = NULL)

\method{spTable}{formula}(dataS, dataP, select, weights = NULL)
}
\arguments{
  \item{dataS}{a \code{data.frame} containing household survey data.}
  \item{dataP}{a \code{data.frame} containing simulated population data.}
  \item{select}{for the formula method, a formula specifying the variables to 
    be used for cross tabulation.  For the default method, an optional 
    character vector defining the columns of \code{dataS} and \code{dataP} to 
    be used.}
  \item{weights}{either a numeric vector containing the (personal) sample 
    weights, or the name of the corresponding column of \code{dataS} (for the 
    default method, the name must be a character string).}
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
  Andreas Alfons
}
\note{
  The class of the argument \code{select} determines the method to be 
  dispatched, not the class of the first argument.
}
\seealso{
  \code{\link{spMosaic}}, \code{\link{tableWt}}
}
\examples{
set.seed(1234)  # for reproducibility
data(eusilcS)   # load sample data
eusilcP <- simStructure(eusilcS)
spTable(eusilcS, eusilcP, select = ~ rb090 + db040, weights = rb050)
}
\keyword{dplot}

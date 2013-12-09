\name{contingencyWt}
\Rdversion{1.1}
\alias{contingencyWt}
\alias{contingencyWt.default}
\alias{contingencyWt.matrix}
\alias{contingencyWt.data.frame}
\title{
  Weighted contingency coefficients
}
\description{
  Compute (weighted) pairwise contingency coefficients.
}
\usage{
contingencyWt(x, \dots)

\S3method{contingencyWt}{default}(x, y, weights = NULL, \dots)

\S3method{contingencyWt}{matrix}(x, weights = NULL, \dots)

\S3method{contingencyWt}{data.frame}(x, weights = NULL, \dots)
}
\arguments{
  \item{x}{for the default method, a vector that can be interpreted as factor.  
    For the matrix and \code{data.frame} methods, the columns should be 
    interpretable as factors.}
  \item{y}{a vector that can be interpreted as factor.}
  \item{weights}{an optional numeric vector containing sample weights.}
  \item{\dots}{for the generic function, arguments to be passed down to the 
    methods, otherwise ignored.}
}
\details{
  The function \code{\link{tableWt}} is used for the computation of the 
  corresponding pairwise contingency tables.
}
\value{
  For the default method, the (weighted) contingency coefficient of \code{x} 
  and \code{y}.

  For the matrix and \code{data.frame} method, a matrix of (weighted) pairwise 
  contingency coefficients for all combinations of columns.  Elements below the 
  diagonal are \code{NA}.
}
\references{
Kendall, M.G. and Stuart, A. (1967) \emph{The Advanced Theory of Statistics, 
Volume 2: Inference and Relationship}. Charles Griffin & Co Ltd, London, 2nd 
edition.}
\author{Andreas Alfons and Stefan Kraft}
\seealso{\code{\link{tableWt}}}
\examples{
data(eusilcS)

## default method
contingencyWt(eusilcS$pl030, eusilcS$pb220a, weights = eusilcS$rb050)

## data.frame method
basic <- c("age", "rb090", "hsize", "pl030", "pb220a")
contingencyWt(eusilcS[, basic], weights = eusilcS$rb050)
}
\keyword{category}

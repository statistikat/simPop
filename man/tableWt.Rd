\name{tableWt}
\Rdversion{1.1}
\alias{tableWt}
\title{
  Weighted cross tabulation
}
\description{
  Compute contingency tables taking into account sample weights.
}
\usage{
tableWt(x, weights = NULL, useNA = c("no", "ifany", "always"))
}
\arguments{
  \item{x}{a vector that can be interpreted as a factor, or a matrix or 
    \code{data.frame} whose columns can be interpreted as factors.}
  \item{weights}{an optional numeric vector containing sample weights.}
  \item{useNA}{a logical indicating whether to include extra \code{NA} levels 
    in the table.}
}
\details{
  For each combination of the variables in \code{x}, the weighted number of 
  occurence is computed as the sum of the corresponding sample weights.  If 
  weights are not specified, the function \code{\link{table}} is applied.
}
\value{
The (weighted) contingency table as an object of class \code{table}, an 
array of integer values.
}
\author{Andreas Alfons and Stefan Kraft}
\seealso{\code{\link{table}}, \code{\link{contingencyWt}}}
\examples{
data(eusilcS)
tableWt(eusilcS[, c("hsize", "db040")], weights = eusilcS$rb050)
tableWt(eusilcS[, c("rb090", "pb220a")], weights = eusilcS$rb050, 
    useNA = "ifany")
}
\keyword{category}

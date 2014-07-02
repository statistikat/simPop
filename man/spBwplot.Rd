\name{spBwplot}
\alias{spBwplot}
\title{
  Weighted box plots
}
\description{
  Produce box-and-whisker plots of continuous or semi-continuous variables, possibly broken down according to conditioning variables and taking into account sample weights.
}
\usage{
spBwplot(inp, x, cond = NULL, horizontal = TRUE, coef = 1.5, zeros = TRUE,
  minRatio = NULL, do.out = FALSE, \dots)
}
\arguments{
  \item{inp}{an object of class \code{linkS4class{synthPopObj}} containing survey sample and synthetic population data.}
  \item{x}{for the default method (currently the only method implemented), a character vector specifying the columns of data available in the sample and the population (specified in input object 'inp') to be plotted.}
  \item{cond}{an optional character vector (of length 1, if used) specifying the conditioning variable.}
  \item{horizontal}{a logical indicating whether the boxes should be horizontal or vertical.}
  \item{coef}{a numeric value that determines the extension of the whiskers.}
  \item{zeros}{a logical indicating whether the variables specified by \code{x} are semi-continuous, i.e., contain a considerable amount of zeros. If \code{TRUE}, the box widths correspond to the proportion of non-zero data points and the (weighted) box plot statistics are computed for these non-zero data points only.}
  \item{minRatio}{a numeric value in \eqn{(0,1]}; if \code{zeros} is \code{TRUE}, the boxes may become unreadable for a large proportion of zeros. In such a case, this can be used to specify a minimum ratio for the box widths. Variable box widths for semi-continuous variables can be suppressed by setting this value to 1.}
  \item{do.out}{a logical indicating whether data points that lie beyond the extremes of the whiskers should be plotted. Note that this is \code{FALSE} by default.}
  \item{\dots}{for the generic function, further arguments to be passed down to methods. For the default method, further arguments to be passed to \code{\link[lattice:xyplot]{bwplot}}.}
}
\details{
 Missing values are ignored for producing box plots and weights are directly extracted from the input object \code{inp}.
}
\value{
  An object of class \code{"trellis"}, as returned by \code{\link[lattice:xyplot]{bwplot}}.
}
\author{
  Andreas Alfons and Bernhard Meindl
}
\seealso{
  \code{\link{spBwplotStats}}, \code{\link[lattice:xyplot]{bwplot}}
}
\examples{
\dontrun{
data(eusilcS)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
synthPop <- simStructure(data=inp, method="direct",
  basicHHvars=c("age", "rb090", "hsize", "pl030", "pb220a"))

# multinomial model with random draws
eusilcM <- simContinuous(synthPop, additional="netIncome", upper=200000, equidist=FALSE)
class(eusilcM)

spBwplot(inp=inp, x="netIncome", cond=NULL)
spBwplot(inp=inp, x="netIncome", cond="rb090", layout=c(1,2))
}
}
\keyword{hplot}

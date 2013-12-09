\name{spBwplot}
\Rdversion{1.1}
\alias{spBwplot}
\alias{spBwplot.default}
\title{
  Weighted box plots
}
\description{  
  Produce box-and-whisker plots of continuous or semi-continuous variables, 
  possibly broken down according to conditioning variables and taking into 
  account sample weights.
}
\usage{
spBwplot(x, \dots)

\method{spBwplot}{default}(x, weights = NULL, cond = NULL, dataS, dataP = NULL, 
         horizontal = TRUE, coef = 1.5, zeros = TRUE, 
         minRatio = NULL, do.out = FALSE, \dots)
}
\arguments{
  \item{x}{for the default method (currently the only method implemented), a 
    character vector specifying the columns of \code{dataS} and \code{dataP} to 
    be plotted.}
  \item{weights}{a character string specifying the column of \code{dataS} that 
    contains the (personal) sample weights.}
  \item{cond}{an optional character vector specifying conditioning variables.}
  \item{dataS}{a \code{data.frame} containing household survey data.}
  \item{dataP}{optional; a \code{data.frame} containing simulated population 
    data or a list of such \code{data.frame}s.}
  \item{horizontal}{a logical indicating whether the boxes should be horizontal 
    or vertical.}
  \item{coef}{a numeric value that determines the extension of the whiskers.}
  \item{zeros}{a logical indicating whether the variables specified by \code{x} 
    are semi-continuous, i.e., contain a considerable amount of zeros.  If 
    \code{TRUE}, the box widths correspond to the proportion of non-zero data 
    points and the (weighted) box plot statistics are computed for these 
    non-zero data points only.}
  \item{minRatio}{a numeric value in \eqn{(0,1]}; if \code{zeros} is 
    \code{TRUE}, the boxes may become unreadable for a large proportion of 
    zeros.  In such a case, this can be used to specify a minimum ratio for the 
    box widths.  Variable box widths for semi-continuous variables can be 
    suppressed by setting this value to 1.}
  \item{do.out}{a logical indicating whether data points that lie beyond the 
    extremes of the whiskers should be plotted.  Note that this is \code{FALSE} 
    by default.}
  \item{\dots}{for the generic function, further arguments to be passed down to 
    methods.  For the default method, further arguments to be passed to 
    \code{\link[lattice:xyplot]{bwplot}}.}
}
\details{
  Missing values are ignored for producing box plots.
}
\value{
  An object of class \code{"trellis"}, as returned by 
  \code{\link[lattice:xyplot]{bwplot}}.
}
\author{
  Andreas Alfons
}
\note{
  A formula interface may be added in the future.
}
\seealso{
  \code{\link{spBwplotStats}}, \code{\link[lattice:xyplot]{bwplot}}
}
\examples{
\dontrun{

## these take some time and are not run automatically
## copy & paste to the R command line

set.seed(1234)  # for reproducibility
data(eusilcS)   # load sample data
# multinomial model with random draws
eusilcM <- simEUSILC(eusilcS, upper = 200000, equidist = FALSE)
# two-step regression
eusilcT <- simEUSILC(eusilcS, method = "twostep")
# plot results by gender
spBwplot("netIncome", "rb050", "rb090", dataS = eusilcS, 
    dataP = list(M = eusilcM, T = eusilcT), layout = c(1, 2))
}
}
\keyword{hplot}

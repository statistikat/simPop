\name{simComponents}
\Rdversion{1.1}
\alias{simComponents}
\title{
  Simulate components of continuous variables of population data
}
\description{
  Simulate components of continuous variables of population data by resampling 
  fractions from survey data.  The continuous variable to be split and any 
  categorical conditioning variables need to be simulated beforehand.
}
\usage{
simComponents(dataS, dataP, w = "rb050", total = "netIncome", 
              components = c("py010n", "py050n", "py090n", 
                "py100n", "py110n", "py120n", "py130n", "py140n"),
              conditional = c(getCatName(total), "pl030"), 
              replaceEmpty = c("sequential", "min"), seed)
}
\arguments{
  \item{dataS}{a \code{data.frame} containing household survey data.}
  \item{dataP}{a \code{data.frame} containing the simulated population data.}
  \item{w}{a character string specifying the column of \code{dataS} that 
    contains the (personal) sample weights, which are used as probability 
    weights for resampling.}
  \item{total}{a character string specifying the continuous variable of 
    \code{dataP} that should be split into components.  Currently, only one 
    variable can be split at a time.}
  \item{components}{a character vector specifying the components in 
    \code{dataS} that should be simulated for the population data.}
  \item{conditional}{an optional character vector specifying categorical 
    contitioning variables for resampling.  The fractions occurring in 
    \code{dataS} are then drawn from the respective subsets defined by these 
    variables.}
  \item{replaceEmpty}{a character string; if \code{conditional} specifies at 
    least two conditioning variables, this determines how replacement cells for 
    empty subsets in the sample are obtained.  If \code{"sequential"}, the 
    conditioning variables are browsed sequentially such that replacement cells 
    have the same value in one conditioning variable and minimum Manhattan 
    distance in the other conditioning variables.  If no such cells exist, 
    replacement cells with minimum overall Manhattan distance are selected.  
    The latter is always done if this is \code{"min"} or only one conditioning 
    variable is used.}
  \item{seed}{optional; an integer value to be used as the seed of the random 
    number generator, or an integer vector containing the state of the random 
    number generator to be restored.}
}
\value{
  A \code{data.frame} containing the simulated population data including the 
  components of the continuous variable specified by \code{total}.
}
\author{Stefan Kraft and Andreas Alfons}
\note{
  The basic household structure, any categorical conditioning variables and the 
  continuous variable to be split need to be simulated beforehand with the 
  functions \code{\link{simStructure}}, \code{\link{simCategorical}} and 
  \code{\link{simContinuous}}.
}
\seealso{
  \code{\link{simStructure}}, \code{\link{simCategorical}}, 
  \code{\link{simContinuous}}, \code{\link{simEUSILC}}
}
\examples{
\dontrun{

## these take some time and are not run automatically
## copy & paste to the R command line

set.seed(1234)  # for reproducibility
data(eusilcS)   # load sample data
eusilcP <- simStructure(eusilcS)
eusilcP <- simCategorical(eusilcS, eusilcP)
basic <- c("age", "rb090", "hsize", "pl030", "pb220a")
eusilcP <- simContinuous(eusilcS, eusilcP, 
    basic = basic, upper = 200000, equidist = FALSE)
# categorize net income for use as conditioning variable
breaks <- getBreaks(eusilcS$netIncome, eusilcS$rb050, 
    upper=Inf, equidist = FALSE)
eusilcS$netIncomeCat <- getCat(eusilcS$netIncome, breaks)
eusilcP$netIncomeCat <- getCat(eusilcP$netIncome, breaks)
# simulate net income components
eusilcP <- simComponents(eusilcS, eusilcP)
summary(eusilcP)
}
}
\keyword{datagen}

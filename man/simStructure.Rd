\name{simStructure}
\Rdversion{1.1}
\alias{simStructure}
\title{
  Simulate the household structure of population data
}
\description{
  Simulate basic categorical variables that define the household structure 
  (typically household ID, age and gender) of population data by resampling 
  from survey data.
}
\usage{
simStructure(dataS, method = c("direct", "multinom", "distribution"), seed)
}
\arguments{
  \item{dataS}{a \code{data.frame} containing household survey data.}
  \item{method}{a character string specifying the method to be used for 
    simulating the household sizes.  Accepted values are \code{"direct"} 
    (estimation of the population totals for each combination of stratum and 
    household size using the Horvitz-Thompson estimator), \code{"multinom"} 
    (estimation of the conditional probabilities within the strata using a 
    multinomial log-linear model and random draws from the resulting 
    distributions), or \code{"distribution"} (random draws from the observed 
    conditional distributions within the strata).}
  \item{seed}{optional; an integer value to be used as the seed of the random 
    number generator, or an integer vector containing the state of the random 
    number generator to be restored.}
}
\value{
  An object of class \code{synthPopObj} containing the simulated population household structure as well as the underlying sample.
}
\author{Bernhard Meindl}
\note{
  The function \code{\link{sample}} is used, which gives results incompatible 
  with those from \R < 2.2.0 and produces a warning the first time this happens 
  in a session.
}
\seealso{
  \code{\link{simCategorical}}, \code{\link{simContinuous}}, 
  \code{\link{simComponents}}, \code{\link{simEUSILC}}
}
\examples{
data(eusilcS)   # load sample data
samp <- specify_sample(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090", additional=c("age","rb090"))
eusilcP <- simStructure(data=samp, method="direct")
summary(eusilcP)
}
\keyword{datagen}

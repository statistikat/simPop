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
simStructure(dataS, hid = "db030", w = "db090", 
             hsize = NULL, strata = "db040", pid = NULL, 
             additional = c("age", "rb090"), 
             method = c("direct", "multinom", "distribution"), 
             keep = TRUE, seed)
}
\arguments{
  \item{dataS}{a \code{data.frame} containing household survey data.}
  \item{hid}{a character string specifying the column of \code{dataS} that 
    contains the household ID.}
  \item{w}{a character string specifying the column of \code{dataS} that 
    contains the (household) sample weights, which are used as probability 
    weights for resampling.}
  \item{hsize}{an optional character string specifying a column of 
    \code{dataS} that contains the household size.  If \code{NULL}, the 
    household sizes are computed.}
  \item{strata}{a character string specifying the column of \code{dataS} that 
    define strata.  Note that this is currently a required argument and only 
    one stratification variable is supported.}
  \item{pid}{an optional character string specifying a column of \code{dataS} 
    that contains the personal ID.}
  \item{additional}{a character vector specifying additional categorical 
    variables of \code{dataS} that define the household structure, typically 
    age and gender.}
  \item{method}{a character string specifying the method to be used for 
    simulating the household sizes.  Accepted values are \code{"direct"} 
    (estimation of the population totals for each combination of stratum and 
    household size using the Horvitz-Thompson estimator), \code{"multinom"} 
    (estimation of the conditional probabilities within the strata using a 
    multinomial log-linear model and random draws from the resulting 
    distributions), or \code{"distribution"} (random draws from the observed 
    conditional distributions within the strata).}
  \item{keep}{a logical indicating whether the original IDs of the 
    corresponding households in the underlying sample should be stored as a 
    variable in the resulting population data.  If \code{TRUE}, the 
    corresponding column name is given by \code{hid} with postfix 
    \code{"Sample"}.}
  \item{seed}{optional; an integer value to be used as the seed of the random 
    number generator, or an integer vector containing the state of the random 
    number generator to be restored.}
}
\value{
  A \code{data.frame} containing the simulated population household structure.
}
\author{Andreas Alfons and Stefan Kraft}
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
set.seed(1234)  # for reproducibility
data(eusilcS)   # load sample data
eusilcP <- simStructure(eusilcS)
summary(eusilcP)
}
\keyword{datagen}

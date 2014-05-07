\name{simCategorical}
\alias{simCategorical}
\title{
  Simulate categorical variables of population data
}
\description{
  Simulate categorical variables of population data.  The household structure 
  of the population data needs to be simulated beforehand.
}
\usage{
simCategorical(synthPopObj, basic, additional, 
  method=c("multinom", "distribution", "ctree", "naivebayes"), 
  limit=NULL, censor=NULL, maxit=500, MaxNWts=1500, eps=NULL, seed=1)
}
\arguments{
  \item{synthPopObj}{a \code{synthPopObj} containing population and household survey data as well as optionally margins in standardized format.}
  \item{basic}{a character vector specifying the columns of \code{dataS} and 
    \code{dataP}, respectively, that define the household structure, typically 
    age, gender and household size.  The default value is \code{c("age", 
    "rb090", "hsize")} if \code{method} is \code{"multinom"}, and 
    \code{c("age", "rb090")} if \code{method} is \code{"distribution"}.}
  \item{additional}{a character vector specifying additional categorical 
    variables of \code{dataS} that should be simulated for the population data.}
  \item{method}{a character string specifying the method to be used for 
    simulating the additional categorical variables.  Accepted values are 
    \code{"multinom"} (estimation of the conditional probabilities using 
    multinomial log-linear models and random draws from the resulting 
    distributions), or \code{"distribution"} (random draws from the observed 
    conditional distributions of their multivariate realizations), \code{ctree} 
    (classification trees) and \code{naivebayes} (classification based on naive 
    bayes approach).}
  \item{limit}{if \code{method} is \code{"multinom"}, this can be used to 
    account for structural zeros.  If only one additional variable is 
    requested, a named list of lists should be supplied.  The names of the list 
    components specify the predictor variables for which to limit the possible 
    outcomes of the response.  For each predictor, a list containing the 
    possible outcomes of the response for each category of the predictor can be 
    supplied.  The probabilities of other outcomes conditional on combinations 
    that contain the specified categories of the supplied predictors are set 
    to 0.  If more than one additional variable is requested, such a list of 
    lists can be supplied for each variable as a component of yet another list, 
    with the component names specifying the respective variables.}
  \item{censor}{if \code{method} is \code{"multinom"}, this can be used to 
    account for structural zeros.  If only one additional variable is 
    requested, a named list of lists or \code{data.frame}s should be 
    supplied.  The names of the list components specify the categories that 
    should be censored.  For each of these categories, a list or 
    \code{data.frame} containing levels of the predictor variables can be 
    supplied.  The probability of the specified categories is set to 0 for the 
    respective predictor levels.  If more than one additional variable is 
    requested, such a list of lists or \code{data.frame}s can be supplied for 
    each variable as a component of yet another list, with the component names 
    specifying the respective variables.}
  \item{maxit, MaxNWts}{control parameters to be passed to 
    \code{\link[nnet]{multinom}} and \code{\link[nnet]{nnet}}.  See the help 
    file for \code{\link[nnet]{nnet}}.}
  \item{eps}{a small positive numeric value, or \code{NULL} (the default).  In 
    the former case and if \code{method} is \code{"multinom"}, estimated 
    probabilities smaller than this are assumed to result from structural zeros 
    and are set to exactly 0.}
  \item{seed}{optional; an integer value to be used as the seed of the random 
    number generator, or an integer vector containing the state of the random 
    number generator to be restored.}
}
\value{
  A \code{data.frame} containing the simulated population data including the 
  categorical variables specified by \code{additional}.
}
\author{Bernhard Meindl and Andreas Alfons and Stefan Kraft}
\note{
  The basic household structure needs to be simulated beforehand with the 
  function \code{\link{simStructure}}.
}
\seealso{
  \code{\link{simStructure}}, \code{\link{simRelation}}, 
  \code{\link{simContinuous}}, \code{\link{simComponents}}
}
\examples{
\dontrun{

## these take some time and are not run automatically
## copy & paste to the R command line

set.seed(1234)  # for reproducibility
data(eusilcS)   # load sample data
eusilcP <- simStructure(eusilcS)
eusilcP <- simCategorical(eusilcS, eusilcP)
summary(eusilcP)
}
}
\keyword{datagen}

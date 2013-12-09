\name{simRelation}
\Rdversion{1.1}
\alias{simRelation}
\title{
  Simulate categorical variables of population data
}
\description{
  Simulate categorical variables of population data taking relationships 
  between household members into account.  The household structure of the 
  population data needs to be simulated beforehand.
}
\usage{
simRelation(dataS, dataP, hid = "hhid", w = "weight", 
            strata = "region", basic = c("age", "sex", "hsize"), 
            relation = "relate", head = "head", direct = NULL, 
            additional = c("nation", "ethnic", "religion"), 
            limit = NULL, censor = NULL, maxit = 500, 
            MaxNWts = 2000, eps = NULL, seed)
}
\arguments{
  \item{dataS}{a \code{data.frame} containing household survey data.}
  \item{dataP}{a \code{data.frame} containing the simulated population 
    household structure.}
  \item{hid}{a character string specifying the columns of \code{dataS} and 
    \code{dataP}, respectively, that contain the household ID.}
  \item{w}{a character string specifying the column of \code{dataS} that 
    contains the (personal) sample weights.}
  \item{strata}{a character string specifying the columns of \code{dataS} and 
    \code{dataP}, respectively, that define strata.  The values are simulated 
    for each stratum separately.  Note that this is currently a required 
    argument and only one stratification variable is supported.}
  \item{basic}{a character vector specifying the columns of \code{dataS} and 
    \code{dataP}, respectively, that define the household structure, typically 
    age, gender and household size.}
  \item{relation}{a character string specifying the columns of \code{dataS} and 
    \code{dataP}, respectively, that define the relationships between the 
    household members.}
  \item{head}{a character string specifying the category of the variable given  
    by \code{relation} that identifies the household head.}
  \item{direct}{a character string specifying categories of the variable given 
    by \code{relation}.  Simulated individuals with those categories directly 
    inherit the values of the additional variables from the household head.  
    The default is \code{NULL} such that no individuals directly inherit values 
    from the household head.}
  \item{additional}{a character vector specifying additional categorical 
    variables of \code{dataS} that should be simulated for the population data.}
  \item{limit}{this can be used to account for structural zeros.  If only 
    one additional variable is requested, a named list of lists should be 
    supplied.  The names of the list components specify the predictor variables 
    for which to limit the possible outcomes of the response.  For each 
    predictor, a list containing the possible outcomes of the response for each 
    category of the predictor can be supplied.  The probabilities of other 
    outcomes conditional on combinations that contain the specified categories 
    of the supplied predictors are set to 0.  If more than one additional 
    variable is requested, such a list of lists can be supplied for each 
    variable as a component of yet another list, with the component names 
    specifying the respective variables.}
  \item{censor}{this can be used to account for structural zeros.  If only 
    one additional variable is requested, a named list of lists or 
    \code{data.frame}s should be supplied.  The names of the list components 
    specify the categories that should be censored.  For each of these 
    categories, a list or \code{data.frame} containing levels of the predictor 
    variables can be supplied.  The probability of the specified categories is 
    set to 0 for the respective predictor levels.  If more than one additional 
    variable is requested, such a list of lists or \code{data.frame}s can be 
    supplied for each variable as a component of yet another list, with the 
    component names specifying the respective variables.}
  \item{maxit, MaxNWts}{control parameters to be passed to 
    \code{\link[nnet]{multinom}} and \code{\link[nnet]{nnet}}.  See the help 
    file for \code{\link[nnet]{nnet}}.}
  \item{eps}{a small positive numeric value, or \code{NULL} (the default).  In 
    the former case, estimated probabilities smaller than this are assumed to 
    result from structural zeros and are set to exactly 0.}
  \item{seed}{optional; an integer value to be used as the seed of the random 
    number generator, or an integer vector containing the state of the random 
    number generator to be restored.}
}
\details{
  The values of a new variable are simulated in three steps, where the second 
  step is optional.  First, the values of the household heads are simulated 
  with multinomial log-linear models.  Second, individuals directly related to 
  the corresponding household head (as specified by the argument \code{direct}) 
  inherit the value of the latter.  Third, the values of the remaining 
  individuals are simulated with multinomial log-linear models in which the 
  value of the respective household head is used as an additional predictor.
}
\value{
  A \code{data.frame} containing the simulated population data including the 
  categorical variables specified by \code{additional}.
}
\author{Andreas Alfons}
\note{
  The basic household structure needs to be simulated beforehand with the 
  function \code{\link{simStructure}}.
}
\seealso{
  \code{\link{simStructure}}, \code{\link{simCategorical}}, 
  \code{\link{simContinuous}}, \code{\link{simComponents}}
}
\examples{
\dontrun{

## these take some time and are not run automatically
## copy & paste to the R command line

set.seed(1234)  # for reproducibility
data(ghanaS)   # load sample data
ghanaP <- simStructure(ghanaS, hid = "hhid", w = "weight", 
    strata = "region", additional = c("age", "sex", "relate"))
ghanaP <- simRelation(ghanaS, ghanaP)
summary(ghanaP)
}
}
\keyword{datagen}

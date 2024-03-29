% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simEUSILC.R
\name{simEUSILC}
\alias{simEUSILC}
\title{Simulate EU-SILC population data}
\usage{
simEUSILC(
  dataS,
  hid = "db030",
  wh = "db090",
  wp = "rb050",
  hsize = NULL,
  strata = "db040",
  pid = NULL,
  age = "age",
  gender = "rb090",
  categorizeAge = TRUE,
  breaksAge = NULL,
  categorical = c("pl030", "pb220a"),
  income = "netIncome",
  method = c("multinom", "twostep"),
  breaks = NULL,
  lower = NULL,
  upper = NULL,
  equidist = TRUE,
  probs = NULL,
  gpd = TRUE,
  threshold = NULL,
  est = "moments",
  const = NULL,
  alpha = 0.01,
  residuals = TRUE,
  components = c("py010n", "py050n", "py090n", "py100n", "py110n", "py120n", "py130n",
    "py140n"),
  conditional = c(getCatName(income), "pl030"),
  keep = TRUE,
  maxit = 500,
  MaxNWts = 1500,
  tol = .Machine$double.eps^0.5,
  nr_cpus = NULL,
  seed
)
}
\arguments{
\item{dataS}{a \code{data.frame} containing EU-SILC survey data.}

\item{hid}{a character string specifying the column of \code{dataS} that
contains the household ID.}

\item{wh}{a character string specifying the column of \code{dataS} that
contains the household sample weights.}

\item{wp}{a character string specifying the column of \code{dataS} that
contains the personal sample weights.}

\item{hsize}{an optional character string specifying a column of
\code{dataS} that contains the household size. If \code{NULL}, the household
sizes are computed.}

\item{strata}{a character string specifying the column of \code{dataS} that
define strata. Note that this is currently a required argument and only one
stratification variable is supported.}

\item{pid}{an optional character string specifying a column of \code{dataS}
that contains the personal ID.}

\item{age}{a character string specifying the column of \code{dataS} that
contains the age of the persons (to be used for setting up the household
structure).}

\item{gender}{a character string specifying the column of \code{dataS} that
contains the gender of the persons (to be used for setting up the household
structure).}

\item{categorizeAge}{a logical indicating whether age categories should be
used for simulating additional categorical and continuous variables to
decrease computation time.}

\item{breaksAge}{numeric; if \code{categorizeAge} is \code{TRUE}, an
optional vector of two or more break points for constructing age categories,
otherwise ignored.}

\item{categorical}{a character vector specifying additional categorical
variables of \code{dataS} that should be simulated for the population data.}

\item{income}{a character string specifying the variable of \code{dataS}
that contains the personal income (to be simulated for the population data).}

\item{method}{a character string specifying the method to be used for
simulating personal income. Accepted values are \code{"multinom"} (for using
multinomial log-linear models combined with random draws from the resulting
ategories) and \code{"twostep"} (for using two-step regression models
combined with random error terms).}

\item{breaks}{if \code{method} is \code{"multinom"}, an optional numeric
vector of two or more break points for categorizing the personal income. If
missing, break points are computed using weighted quantiles.}

\item{lower, upper}{numeric values; if \code{method} is \code{"multinom"} and
\code{breaks} is \code{NULL}, these can be used to specify lower and upper
bounds other than minimum and maximum, respectively. Note that if \code{gpd}
is \code{TRUE} (see below), \code{upper} defaults to \code{Inf}.}

\item{equidist}{logical; if \code{method} is \code{"multinom"} and
\code{breaks} is \code{NULL}, this indicates whether the (positive) default
break points should be equidistant or whether there should be refinements in
the lower and upper tail (see \code{\link{getBreaks}}).}

\item{probs}{numeric vector with values in \eqn{[0, 1]}; if \code{method} is
\code{"multinom"} and \code{breaks} is \code{NULL}, this gives probabilities
for quantiles to be used as (positive) break points. If supplied, this is
preferred over \code{equidist}.}

\item{gpd}{logical; if \code{method} is \code{"multinom"}, this indicates
whether the upper tail of the personal income should be simulated by random
draws from a (truncated) generalized Pareto distribution rather than a
uniform distribution.}

\item{threshold}{a numeric value; if \code{method} is \code{"multinom"},
values for categories above \code{threshold} are drawn from a (truncated)
generalized Pareto distribution.}

\item{est}{a character string; if \code{method} is \code{"multinom"}, the
estimator to be used to fit the generalized Pareto distribution.}

\item{const}{numeric; if \code{method} is \code{"twostep"}, this gives a
constant to be added before log transformation.}

\item{alpha}{numeric; if \code{method} is \code{"twostep"}, this gives
trimming parameters for the sample data. Trimming is thereby done with
respect to the variable specified by \code{additional}. If a numeric vector
of length two is supplied, the first element gives the trimming proportion
for the lower part and the second element the trimming proportion for the
upper part. If a single numeric is supplied, it is used for both. With
\code{NULL}, trimming is suppressed.}

\item{residuals}{logical; if \code{method} is \code{"twostep"}, this
indicates whether the random error terms should be obtained by draws from
the residuals. If \code{FALSE}, they are drawn from a normal distribution
(median and MAD of the residuals are used as parameters).}

\item{components}{a character vector specifying the income components in
\code{dataS} (to be simulated for the population data).}

\item{conditional}{an optional character vector specifying categorical
contitioning variables for resampling of the income components. The
fractions occurring in \code{dataS} are then drawn from the respective
subsets defined by these variables.}

\item{keep}{a logical indicating whether variables computed internally in
the procedure (such as the original IDs of the corresponding households in
the underlying sample, age categories or income categories) should be stored
in the resulting population data.}

\item{maxit, MaxNWts}{control parameters to be passed to
\code{\link[nnet]{multinom}} and \code{\link[nnet]{nnet}}. See the help file
for \code{\link[nnet]{nnet}}.}

\item{tol}{if \code{method} is \code{"twostep"}, a small positive numeric
value or \code{NULL} (see \code{\link{simContinuous}}).}

\item{nr_cpus}{if specified, an integer number defining the number of cpus
that should be used for parallel processing.}

\item{seed}{optional; an integer value to be used as the seed of the random
number generator, or an integer vector containing the state of the random
number generator to be restored.}
}
\value{
An object of class \code{\linkS4class{simPopObj}} containing the
simulated EU-SILC population data as well as the underlying sample.
}
\description{
Simulate population data for the European Statistics on Income and Living
Conditions (EU-SILC).
}
\note{
This is a wrapper calling \code{\link{simStructure}},
\code{\link{simCategorical}}, \code{\link{simContinuous}} and
\code{\link{simComponents}}.
}
\examples{

data(eusilcS) # load sample data

\dontrun{
## long computation time
# multinomial model with random draws
eusilcM <- simEUSILC(eusilcS, upper = 200000, equidist = FALSE
, nr_cpus = 1)
summary(eusilcM)

# two-step regression
eusilcT <- simEUSILC(eusilcS, method = "twostep", nr_cpus = 1)
summary(eusilcT)
}

}
\seealso{
\code{\link{simStructure}}, \code{\link{simCategorical}},
\code{\link{simContinuous}}, \code{\link{simComponents}}
}
\author{
Andreas Alfons and Stefan Kraft and Bernhard Meindl
}
\keyword{datagen}

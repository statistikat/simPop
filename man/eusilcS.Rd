\name{eusilcS}
\alias{eusilcS}
\docType{data}
\title{
Synthetic EU-SILC survey data
}
\description{
This data set is synthetically generated from real Austrian EU-SILC (European 
Union Statistics on Income and Living Conditions) data.
}
\usage{data(eusilcS)}
\format{
  A data frame with 11725 observations on the following 18 variables.
  \describe{
    \item{\code{db030}}{integer; the household ID.}
    \item{\code{hsize}}{integer; the number of persons in the household.}
    \item{\code{db040}}{factor; the federal state in which the household is 
      located (levels \code{Burgenland}, \code{Carinthia}, 
      \code{Lower Austria}, \code{Salzburg}, \code{Styria}, \code{Tyrol}, 
      \code{Upper Austria}, \code{Vienna} and \code{Vorarlberg}).}
    \item{\code{age}}{integer; the person's age.}
    \item{\code{rb090}}{factor; the person's gender (levels \code{male} and 
      \code{female}).}
    \item{\code{pl030}}{factor; the person's economic status (levels 
      \code{1} = working full time, \code{2} = working part time, \code{3} = 
      unemployed, \code{4} = pupil, student, further training or unpaid work 
      experience or in compulsory military or community service, \code{5} = in 
      retirement or early retirement or has given up business, \code{6} = 
      permanently disabled or/and unfit to work or other inactive person, 
      \code{7} = fulfilling domestic tasks and care responsibilities).}
    \item{\code{pb220a}}{factor; the person's citizenship (levels \code{AT}, 
      \code{EU} and \code{Other}).}
    \item{\code{netIncome}}{numeric; the personal net income.}
    \item{\code{py010n}}{numeric; employee cash or near cash income (net).}
    \item{\code{py050n}}{numeric; cash benefits or losses from self-employment 
      (net).}
    \item{\code{py090n}}{numeric; unemployment benefits (net).}
    \item{\code{py100n}}{numeric; old-age benefits (net).}
    \item{\code{py110n}}{numeric; survivor's benefits (net).}
    \item{\code{py120n}}{numeric; sickness benefits (net).}
    \item{\code{py130n}}{numeric; disability benefits (net).}
    \item{\code{py140n}}{numeric; education-related allowances (net).}
    \item{\code{db090}}{numeric; the household sample weights.}
    \item{\code{rb050}}{numeric; the personal sample weights.}
  }
}
\details{
The data set consists of 4641 households and is used as sample data in some of 
the examples in package \code{simPopulation}.  Note that it is included for 
illustrative purposes only.  The sample weights do not reflect the true 
population sizes of Austria and its regions.  The resulting population data is 
about 100 times smaller than the real population size to save computation time.

Only a few of the large number of variables in the original survey are included 
in this example data set.  The variable names are rather cryptic codes, but 
these are the standardized names used by the statistical agencies.  Furthermore, 
the variables \code{hsize}, \code{age} and \code{netIncome} are not included in 
the standardized format of EU-SILC data, but have been derived from other 
variables for convenience.
}
\source{
This is a synthetic data set based on Austrian EU-SILC data from 2006.  The 
original sample was provided by Statistics Austria.
}
\references{
Eurostat (2004) Description of target variables: Cross-sectional and 
longitudinal. \emph{EU-SILC 065/04}, Eurostat.
}
\examples{
data(eusilcS)
summary(eusilcS)
}
\keyword{datasets}

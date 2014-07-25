\name{calibPop}
\alias{calibPop}
\docType{data}
\title{
Calibration of 0/1 weights by Simulated Annealing
}
\description{
A Simulated Annealing Algorithm for calibration of synthetic population data available in a \code{\linkS4class{synthPopObj}}-object. The aims is to find, given a population, a combination of different households which optimally satisfy, in the sense of an acceptable error, a given table of specific known marginals. The known marginals are also already available in slot 'table' of the input object 'inp'.
}
\usage{
calibPop(inp, split, temp=1, eps.factor=0.05, maxiter=200, temp.cooldown=0.9,
  factor.cooldown=0.85, min.temp=10^-3, verbose=FALSE)
}
\arguments{
  \item{inp}{
an object of class \code{\linkS4class{synthPopObj}} with slot 'table' being non-null! (see \code{\link{addKnownMargins}}.
}
  \item{split}{
given strata in which the problem will be split. Has to correspond to a column population data (slot 'pop' of input argument 'inp') . For example \code{split = c("region")}, problem will be split for different regions. Parallel computing is performed automatically, if possible.
}
  \item{temp}{
starting temperatur for simulated annealing algorithm
}
  \item{eps.factor}{
a factor (between 0 and 1) specifying the acceptance error. For example eps.factor = 0.05 results in an acceptance error for the objective function of 0.05*sum(totals)
}
  \item{maxiter}{
maximum iterations during a temperature step.
}
  \item{temp.cooldown}{
a factor (between 0 and 1) specifying the rate at which temperature will be reduced in each step.
}
  \item{factor.cooldown}{
a factor (between 0 and 1) specifying the rate at which the number of permutations of housholds, in each iteration, will be reduced in each step.
}
  \item{min.temp}{
minimal temperature at which the algorithm will stop.
}
  \item{verbose}{
boolean variable; if TRUE some additional verbose output is provided, however only if \code{split} is NULL. Otherwise the computation is performed in parallel and no useful output can be provided.
}
}
\details{
Calibrates data using simulated annealing. The algorithm searches for a (near) optimal combination of different households, by swaping housholds at random in each iteration of each temperature level.
During the algorithm as well as for the output the optimal (or so far best) combination will be indicated by a logical vector containg only 0s (not inculded) and 1s (included in optimal selection).
The objective function for simulated annealing is defined by the sum of absolute differences between target marginals and synthetic marginals (=marginals of synthetic dataset). The sum of target marginals can at most be as large as the sum of target marginals.
For every factor-level in \dQuote{split}, data must at least contain as many entries of this kind as target marginals.

Possible donors are automatically generated within the procedure.
}
\value{
Returns an object of class \code{\linkS4class{synthPopObj}} with an updated population listed in slot 'pop'.
}
\author{
Bernhard Meindl, Johannes Gussenbauer and Matthias Templ
}
\examples{
\dontrun{
data(eusilcS) # load sample data
data(eusilcP) # population data
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
synthPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
synthPop <- simCategorical(synthPop, additional=c("pl030", "pb220a"), method="multinom")

# add margins
margins <- as.data.frame(
  xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$gender + eusilcP$citizenship))
colnames(margins) <- c("db040", "rb090", "pb220a", "freq")
synthPop <- addKnownMargins(synthPop, margins)

# apply simulated annealing
synthPop_adj <- calibPop(synthPop, split="db040", temp=1, eps.factor=0.1)
}
}
\keyword{datasets}

\name{simannealing}
\alias{simannealing}
\docType{data}
\title{
Simulated Annealing to calibrate to known marginals
}
\description{
A Simulated Annealing Algorithm for calibrartion of synthetic microdata. The aims is to find, given a population, a combination of different households which optimaly satisfy, in the sense of an acceptable error, a given table of specific attributes.  
}
\usage{
simannealing(data, totals, hid, parameter, split, temp = 30, eps.factor = 0.05, maxiter=200, temp.cooldown = 0.975, factor.cooldown = 0.85, min.temp = 10^-3, sample = FALSE)
}
\arguments{
  \item{data}{
a data.table object containing microdata. 
%Rows must correspond with different attributes in totals.
}
\item{totals}{
a table object containing the target-marginals. Entries need to display aggregated information on individuals (e.g. age, gender, household size, ect...). Attributes can freely be chosen, but need to correspond to a column in the data. 
}
\item{hid}{
spezifying household-ID column in data
}
\item{parameter}{
character vector containg the variable names to be optimized for.
}
  \item{split}{
given strata in which the problem will be split. Has to correspond to a column in data. For example \code{split = c("region")}, problem will be split for different regions. Parallel computing is automatically applied when working under Linux.
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
  \item{sample}{
boolean variable; if TRUE specific weights will be calculated for a targeted household swapping.
}
}
\details{
Calibrates microdata using simulated annealing. Finds optimal combination of different households, by swaping housholds at random (or with sample weights with sample=TRUE) in each iteration of each temperature level.
During the algorithm as well as for the output the optimal (or so far best) combination will be indicated by a logical vector containg only 0s (not inculded) and 1s (included in optimal selection).
Objective function for simulated annealing is defined by the sum of absolute differences between target marginals and synthetic marginals (=marginals of synthetic microdataset). 
The sum of target marginals can at most be as large as the sum of target marginals. 
For every factor-level in \dQuote{split}, data must at least contain as many entries of this kind as target marginals.
}
\value{
Returns a list including the solution, value of objective function, temperature level, ect...
}
\author{
Bernhard Meindl, Johannes Gussenbauer and Matthias Templ
}
\examples{
\dontrun{
library(data.table)
library(simFrame)
library(simPopulation)
library(parallel)

data(eusilcP)
data(eusilcS)

## synth pop:
pop <- eusilcP
colnames(pop)[3] <- "hhsize"
#eusilcP <- data.table(eusilcP, keep.rownames=TRUE)
## donor data:
donors <- sampHH(pop, strata="region", 
          hsize="hhsize")

## microdata (1) and donor (0)
pop <- rbind(pop, donors)
pop$weight <- rep(c(1,0), c(nrow(eusilcP), nrow(donors)))
pop$new.weights <- pop$weights

index <- pop$weights == 1
tab <- table(pop$region[index],pop$gender[index], pop$hsize[index])

## create target marginals
totals <- tableWt(eusilcS[,c("db040", "rb090","hsize")], weights=eusilcS$rb050)
totals <- nrow(pop[index])/sum(totals) * totals
totals <- ceiling(totals)
totals <- as.data.table(totals)
setnames(totals, c("region","gender","hsize","N"))

data <- data.table(pop, keep.rownames=TRUE)
totals <- totals
hid <- "hid"
parameter <- c("gender","hsize")
split <- "region"
temp <- 10
eps.factor <-  0.05
maxiter <- 250
temp.cooldown <- 0.90
factor.cooldown <- 0.95
min.temp <- 10^-2
sample <- TRUE

system.time(testrun <- simannealing(data, totals, hid, parameter, split, temp, eps.factor, maxiter, temp.cooldown, factor.cooldown, min.temp, sample))
}
}
\keyword{datasets}

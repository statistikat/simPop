\name{simInitSpatial}
\alias{simInitSpatial}
\title{Generation of smaller regions given an existing spatial variable and a table}.
\description{
  This function allows to manipulate an object of class \code{\linkS4class{synthPopObj}} in a way that a new variable containing smaller regions within an already existing broader region is generated. The distribution of the smaller region within the broader region is respected.

The distributional information must be contained in an input table that holds combinations of characteristics of the broader region and the smaller regions as well as population counts (which may be available from a census).
}
\usage{
simInitSpatial(synthPopObj, additional, region, tspatial) 
}
\arguments{
  \item{synthPopObj}{an object of class \code{\linkS4class{synthPopObj}}.}
  \item{additional}{a character vector of length one holding the variable name of the variable containing smaller geographical units. This variable name must be available as a column in input argument \code{tspatial}.}
  \item{region}{a character vector of length one holding the variable name of the broader region. This variable must be available in the input \code{tspatial} as well as in the sample and population slots of input \code{synthPopObj}.}
  \item{tspatial}{a data.frame containing three columns. The broader region (with the variable name being the same as in input \code{region}, the smaller geographical units (with the variable name being the same as in input \code{additional} and a third column containing a numeric vector holding counts.))}
}

\value{
  An object of class \code{\linkS4class{synthPopObj}} with an additional variable in the synthetic population slot.
}
\author{
  Bernhard Meindl
}

\examples{
\dontrun{
data(eusilcS)
data(eusilcP)

# no districts are available in the population, so we have to generate those
# we randomly assign districts within "region" in the eusilc population data
# each hh has the same district
simulate_districts <- function(inp) {
  hhid <- "hid"
  region <- "region"
  
  a <- inp[!duplicated(inp[,hhid]),c(hhid, region)]
  spl <- split(a, a[,region])
  regions <- unique(inp[,region])
  
  tmpres <- lapply(1:length(spl), function(x) {
    codes <- paste(x, 1:sample(3:9,1), sep="")
    spl[[x]]$district <- sample(codes, nrow(spl[[x]]), replace=TRUE)
    spl[[x]]
  })
  tmpres <- do.call("rbind", tmpres)
  tmpres <- tmpres[,-c(2)]
  out <- merge(inp, tmpres, by.x=c(hhid), by.y=hhid, all.x=TRUE)
  invisible(out)
}

eusilcP <- simulate_districts(eusilcP)
table(eusilcP$district)

# we generate a synthetic population
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
synthPopObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))

# we generate the input table using the broad region (variable 'region')
# and the districts, we have generated before.
# we
tab <- as.data.frame(xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$district))
colnames(tab) <- c("db040", "district", "Freq")

synthPopObj <- simInitSpatial(synthPopObj, additional="district", region="db040", tspatial=tab)
}
}
\keyword{manip}

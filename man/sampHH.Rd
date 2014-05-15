\name{sampHH}
\alias{sampHH}
\title{
Sample households from given microdata.
}
\description{
The function samples households from microdata containing personal and 
household information.
}
\usage{
sampHH(pop, sizefactor = 1, hid = "hid", strata = "region", hsize = NULL)
}
\arguments{
\item{pop}{
data frame containing households and persons
}
\item{sizefactor}{
factor of how many times the initial population should be resampled
}
\item{hid}{
string specifying the name of the household-id variable in the data.
}
\item{strata}{
can be used to sample within strata.
}
\item{hsize}{
string specifying the name of the household size variable in the data.
}
}
\details{
households are drawn from the data and new ID's are generated for 
the new households.
}
\value{
the data frame of new households.
}
\author{
Bernhard Meindl and Matthias Templ
}
\examples{
data(eusilcP)
pop <- eusilcP
colnames(pop)[3] <- "hhsize"

system.time(x1 <- sampHH(pop, strata="region", hsize="hhsize"))
dim(x1)
system.time(x1 <- sampHH(pop, sizefactor=4, strata="region", hsize="hhsize"))
dim(x1)
system.time(x2 <- sampHH(pop, strata=NULL, hsize="hhsize"))

pop <- pop[,-which(colnames(pop)=="hhsize")]
system.time(y1 <- sampHH(pop, strata="region", hsize=NULL))
system.time(y2 <- sampHH(pop, strata=NULL, hsize=NULL))
}
\keyword{ manip }

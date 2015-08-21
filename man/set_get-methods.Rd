\name{get_set-methods}
\docType{methods}
\alias{pop}
\alias{pop<-}
\alias{samp}
\alias{samp<-}
\alias{samp,simPopObj-method}
\alias{pop,simPopObj-method}
\alias{samp<-,simPopObj-method}
\alias{pop<-,simPopObj-method}

\title{Extract and modify variables from population or sample data stored in an object of class \code{\link{simPopObj-class}}.}
\usage{
samp(obj, var)

samp(obj, var) <- value

pop(obj, var)

pop(obj, var) <- value
}
\arguments{
\item{obj}{An object of class \code{\link{simPopObj-class}}}

\item{var}{variable name or index for the variable in slot 'samp' of object with the slot name to be
accessed.}

\item{value}{Content replacing whatever the variable in slot \code{var} in \code{obj}
currently holds.}
}
\value{
Returns an object of class \code{\link{simPopObj-class}} with the
appropriate replacement.
}
\description{
Using \code{\link{samp}} \code{\link{samp<-}} it is possible to extract or rather modify variables
of the sample data within slot \code{data} in slot \code{sample} of the \code{\link{simPopObj-class}}-object.
Using \code{\link{pop}} \code{\link{pop<-}} it is possible to extract or rather modify variables
of the synthetic population within in slot \code{data} in slot \code{sample} of the \code{\link{simPopObj-class}}-object.
}
\examples{
data(eusilcS)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040",
weight="db090")
simPopObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))

## get/set variables in sample-object of simPopObj
head(samp(simPopObj, var="age"))
samp(simPopObj, var="newVar") <- 1
head(samp(simPopObj, var="newVar"))
## deleting is also possible
samp(simPopObj, var="newvar") <- NULL
head(samp(simPopObj, var="newvar"))
## extract multiple variables
head(samp(simPopObj, var=c("db030","db040")))

## get/set variables in pop-object of simPopObj
head(pop(simPopObj, var="age"))
pop(simPopObj, var="newVar") <- 1
head(pop(simPopObj, var="newVar"))
## deleting is also possible
pop(simPopObj, var="newvar") <- NULL
head(pop(simPopObj, var="newvar"))
## extract multiple variables
head(pop(simPopObj, var=c("db030","db040")))
}
\author{
Bernhard Meindl
}
\seealso{
\code{\link{simPopObj-class}},\code{\link{pop}}, \code{\link{pop<-}}, \code{\link{samp<-}}, \code{\link{manageSimPopObj}}
}
\keyword{manip}
\keyword{methods}


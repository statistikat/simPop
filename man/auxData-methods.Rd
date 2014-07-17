\name{auxData-methods}
\docType{methods}
\alias{auxData}
\alias{auxData-methods}
\alias{auxData<-}
\alias{auxData<--methods}
\alias{auxData,synthPopObj-method}
\alias{auxData<-,synthPopObj,dataObj-method}
\title{Methods for Function \code{auxData}}
\description{
 allows to query or replace slot 'sample' of an \code{\linkS4class{synthPopObj}} with an object of class \code{\linkS4class{dataObj}}.
}
\section{Methods}{
\describe{
\item{\code{signature(object="synthPopObj")}}{
\code{auxData(inp)} returns slot 'sample' of \code{\linkS4class{synthPopObj}}-object 'iinp'.
}
\item{\code{signature(object="synthPopObj", value="dataObj")}}{
\code{auxData(inp) <- x} replaces slot 'sample' of input \code{inp} with \code{x}. \code{x} must be an object of class \code{\linkS4class{dataObj}}.
}}}

\examples{
data(eusilcS)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
inpObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
class(inpObj)
str(auxData(inpObj))

inp2 <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="pb220a", weight="db090")
auxData(inpObj) <- inp2
str(auxData(inpObj))
}
\keyword{methods}

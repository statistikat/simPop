\name{simPopMethods}
\docType{methods}
\alias{sampleObj}
\alias{sampleObj-methods}
\alias{sampleObj<-}
\alias{sampleObj<--methods}
\alias{sampleData}
\alias{popObj}
\alias{popObj-methods}
\alias{popObj<-}
\alias{popObj<--methods}
\alias{popData}
\alias{tableObj}
\alias{tableObj-methods}
\alias{sampleObj,simPopObj-method}
\alias{sampleObj<-,simPopObj,dataObj-method}
\alias{sampleData,simPopObj-method}
\alias{popObj,simPopObj-method}
\alias{popObj<-,simPopObj,dataObj-method}
\alias{popData,simPopObj-method}
\alias{tableObj,simPopObj-method}

\title{Methods to query/set slots from objects of class \code{\linkS4class{simPopObj}}.}
\description{
 The methods allow to query or replace various information of \code{\linkS4class{simPopObj}}-objects.
}
\section{Methods}{
\describe{
  Below, the methods to access/set various slots of \code{\linkS4class{simPopObj}}-objects are shown.
\item{\code{sampleObj} with \code{signature(object="simPopObj")}}{
  returns slot 'sample' of the input object.
}
\item{\code{sampleObj<-} with \code{signature(object="simPopObj", value="dataObj")}}{
  replaces slot 'sample' of the input object with the object \code{value}.
}
\item{\code{sampleData} with \code{signature(object="simPopObj")}}{
  returns slot 'data' of slot 'sample' of the input object.
}
\item{\code{popObj} with \code{signature(object="simPopObj")}}{
  returns slot 'pop' of the input object.
}
\item{\code{popObj<-} with \code{signature(object="simPopObj", value="dataObj")}}{
  replaces slot 'pop' of the input object with the object \code{value}.
}
\item{\code{popData} with \code{signature(object="simPopObj")}}{
  returns slot 'data' of slot 'pop' of the input object.
}
\item{\code{tableObj} with \code{signature(object="simPopObj")}}{
  returns slot 'table' of the input object.
}
}
}

\examples{
data(eusilcS)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
inpObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
class(inpObj)
str(sampleObj(inpObj))

inp2 <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="pb220a", weight="db090")
sampleObj(inpObj) <- inp2
str(sampleObj(inpObj))
}
\keyword{methods}

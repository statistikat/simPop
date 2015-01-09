\name{addWeights-methods}
\docType{methods}
\alias{addWeights-methods}
\alias{addWeights<-}
\alias{addWeights<--methods}
\alias{addWeights<-,dataObj,list-method}
\alias{addWeights<-,simPopObj,list-method}
\title{Methods for Function \code{addWeights}}
\description{
 allows to modify sampling weights of an \code{\linkS4class{dataObj}} or \code{\linkS4class{simPopObj}}-object.
As input the output of \code{\link{calibSample}} must be used.
}
\section{Methods}{
\describe{
\item{\code{signature(object="dataObj", value="list")}}{
\code{addWeight(inp) <- x} replaces the variable in slot "data" containing sampling weights (specified by slot "weight" of input \code{inp} with \code{x}. \code{x} must be a list that was computed using \code{\link{calibWeights}}.
}
\item{\code{signature(object="simPopObj", value="list")}}{
\code{addWeight(inp) <- x} updates slot "sample" of the provided object of class \code{\linkS4class{simPopObj}}. It replaces the variable in slot "data" that contains sampling weights (specified by slot "weight") of this slot with \code{x}. \code{x} must be a list that was computed using \code{\link{calibWeights}}.
}}}

\examples{
\dontrun{
data(eusilcS)
data(totalsRG)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
addWeights(inp) <- calibSample(inp, totalsRG)
}
}
\keyword{methods}

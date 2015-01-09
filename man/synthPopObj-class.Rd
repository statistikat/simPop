\name{simPopObj-class}
\docType{class}
\alias{simPopObj-class}
\alias{show,simPopObj-method}

\title{Class \code{"simPopObj"}}
\description{
An object that is used throughout the package containing information on the sample (in slot \code{sample}), the population (slot \code{pop}) and optionally some margins in form of a table (slot \code{table}).
}
\section{Objects from the Class}{
Objects are automatically created in function \code{\link{simStructure}}.
}
\section{Slots}{
  \describe{
    \item{\code{sample}:}{Object of class \code{dataObjOrNULL} containing sample information.}
    \item{\code{table}:}{Object of class \code{dataframeOrNULL} containing population margins.}
    \item{\code{pop}:}{Object of class \code{dataObjOrNULL} containing a (simulated) population.}
    \item{\code{basicHHvars}:}{character vector specifying variables available in both \code{sample} and \code{population} that define the (basic) household structure.}
  }
}
\section{Methods}{
A show method is defined with class "simPopObj".
}
\author{
  Bernhard Meindl and Matthias Templ
}
\seealso{
\code{\linkS4class{dataObj}}
}
\examples{
showClass("simPopObj")

## show method: generate an object of class simPop first
data(eusilcS)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
eusilcP <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
class(eusilcP)
## shows some basic information:
eusilcP
}
\keyword{classes}

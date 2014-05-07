\name{synthPopObj-class}
\docType{class}
\alias{synthPopObj-class}

\title{Class \code{"synthPopObj"}}
\description{
An object that is used throughout the package containing information on the sample (in slot \code{sample}), the population (slot \code{pop}) and optionally some margins in form of a table (slot \code{table}).
}
\section{Objects from the Class}{
Objects are automatically created in function \code{\link{simStructure}}.
}
\section{Slots}{
  \describe{
    \item{\code{sample}:}{Object of class \code{"sampleObjOrNULL"} containing sample information }
    \item{\code{table}:}{Object of class \code{"dataframeOrNULL"} containing population margins. }
    \item{\code{pop}:}{Object of class \code{"popObjOrNULL"} containing a (simulated) population }
  }
}
\section{Methods}{
No methods defined with class "synthPopObj" in the signature.
}
\author{
  Bernhard Meindl
}
\seealso{
\code{\linkS4class{sampleObj}} \code{\linkS4class{popObj}}
}
\examples{
showClass("synthPopObj")
}
\keyword{classes}

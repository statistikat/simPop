\name{popObj-class}
\Rdversion{1.1}
\docType{class}
\alias{popObj-class}

\title{Class \code{"popObj"}}
\description{
Objects of this class contain information on a population, similarily to objects of class \code{\linkS4class{sampleObj}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("popObj", ...)} but are usually automatically created when using \code{\link{simStructure}}.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"dataframeOrNULL"} containing the survey data.}
    \item{\code{hhid}:}{Object of class \code{"characterOrNULL"} containing variable name holding household ids in \code{data}. }
    \item{\code{hhsize}:}{Object of class \code{"characterOrNULL"} containing variable name holding household sizes in \code{data}. }
    \item{\code{pid}:}{Object of class \code{"characterOrNULL"} containing variable name holding person ids in \code{data}. }
    \item{\code{strata}:}{Object of class \code{"characterOrNULL"} containing variable name holding strata information in \code{data}. }
    \item{\code{additional}:}{Object of class \code{"characterOrNULL"} containing variable names holding additional important information. }
  }
}
\section{Methods}{
No methods defined with class "popObj" in the signature.
}
\author{
Bernhard Meindl
}

\seealso{
\code{\linkS4class{synthPopObj}} \code{\linkS4class{popObj}}
}
\examples{
showClass("popObj")
}
\keyword{classes}

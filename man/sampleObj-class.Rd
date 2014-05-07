\name{sampleObj-class}
\docType{class}
\alias{sampleObj-class}

\title{Class \code{"sampleObj"}}
\description{
This class contains several useful information that is receycled throughout using functionality from this package. The household-based survey data itself are stored in slot \code{data} while several variable names containing useful information on weights, household- and person-ids, household sizes are stored in respective slots.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sampleObj", ...)} or by using \code{\link{specify_sample}}.
}
\section{Slots}{
  \describe{
    \item{\code{weight}:}{Object of class \code{"characterOrNULL"} containing variable name holding sampling weights in \code{data}.}
    \item{\code{data}:}{Object of class \code{"dataframeOrNULL"} containing the survey data.}
    \item{\code{hhid}:}{Object of class \code{"characterOrNULL"} containing variable name holding household ids in \code{data}. }
    \item{\code{hhsize}:}{Object of class \code{"characterOrNULL"} containing variable name holding household sizes in \code{data}. }
    \item{\code{pid}:}{Object of class \code{"characterOrNULL"} containing variable name holding person ids in \code{data}. }
    \item{\code{strata}:}{Object of class \code{"characterOrNULL"} containing variable name holding strata information in \code{data}. }
    \item{\code{additional}:}{Object of class \code{"characterOrNULL"} containing variable names holding additional important information. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{popObj}"}, directly.
}
\section{Methods}{
No methods defined with class "sampleObj" in the signature.
}
\author{
Bernhard Meindl
}

\seealso{
\code{\linkS4class{synthPopObj}} \code{\linkS4class{popObj}}
}
\examples{
showClass("sampleObj")
}
\keyword{classes}

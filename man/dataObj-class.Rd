\name{dataObj-class}
\docType{class}
\alias{dataObj-class}

\title{Class \code{"dataObj"}}
\description{
Objects of this class contain information on a population or survey.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("dataObj", ...)} but are usually automatically created when using \code{\link{simStructure}}.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{dataframeOrNULL} containing the survey data.}
    \item{\code{hhid}:}{Object of class \code{characterOrNULL} containing variable name holding household ids in \code{data}. }
    \item{\code{pid}:}{Object of class \code{characterOrNULL} containing variable name holding person ids in \code{data}. }
    \item{\code{hhsize}:}{Object of class \code{characterOrNULL} containing variable name holding household sizes in \code{data}. }
    \item{\code{weight}:}{Object of class \code{characterOrNULL} containing weights in \code{data}. }
    \item{\code{strata}:}{Object of class \code{characterOrNULL} containing variable name holding strata information in \code{data}. }
    \item{\code{ispopulation}:}{Object of class \code{logicalOrNULL} defining if slot \code{data} contains population data (TRUE) or sample data (FALSE).}
  }
}
\section{Methods}{
No methods defined with class "dataObj" in the signature.
}
\author{
Bernhard Meindl
}

\seealso{
\code{\linkS4class{synthPopObj}}
}
\examples{
showClass("dataObj")
}
\keyword{classes}

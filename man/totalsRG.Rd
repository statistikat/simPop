\name{totalsRG}
\alias{totalsRG}
\docType{data}
\title{
Population totals Region times Gender for Austria 2006
}
\description{
Population characteristics Region times Gender from Austria.
}
\usage{data(popRegionGender)}
\format{
  totalsRG: A data frame with 18 observations on the following 3 variables.
  \describe{
    \item{\code{rb090}}{gender; a factor with levels \code{female} \code{male}}
    \item{\code{db040}}{region; a factor with levels \code{Burgenland} \code{Carinthia} \code{Lower Austria,} \code{Salzburg} \code{Styria} \code{Tyrol} \code{Upper Austria} \code{Vienna} \code{Vorarlberg}}
    \item{\code{Freq}}{totals; a numeric vector}
  }
  totalsRGtab: a two-dimensional table holding the same information
}
\source{
StatCube - statistical data base, \url{http://www.statistik.at/web_de/services/datenbank_superstar/}
}
\examples{
data(totalsRG)
totalsRG
data(totalsRGtab)
totalsRGtab
}
\keyword{datasets}

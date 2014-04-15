\name{ipu}
\alias{ipu}
\title{
  iterative proportional updating
}
\description{
  adjust sampling weights to given totals based on household-level and/or individual level constraints
}
\usage{
ipu(inp, con, hid=NULL, eps=1e-07, verbose=FALSE)
}
\arguments{
  \item{inp}{a \code{data.frame} or \code{data.table} containing holding counts for household and/or personal level attributes that should be fitted.}
  \item{con}{named list with each list element holding a constraint total with list-names relating to column-names in \code{inp}.}
  \item{hid}{character vector specifying the variable containing household-ids within \code{inp} or NULL if such a variable does not exist.}
  \item{eps}{number specifiying convergence limit} 
  \item{verbose}{if TRUE, ipu will print some progress information.} 
}

\examples{
inp <- as.data.frame(matrix(0, nrow=8, ncol=6))
colnames(inp) <- c("hhid","hh1","hh2","p1","p2","p3")
inp$hhid <- 1:8
inp$hh1[1:3] <- 1
inp$hh2[4:8] <- 1
inp$p1 <- c(1,1,2,1,0,1,2,1)
inp$p2 <- c(1,0,1,0,2,1,1,1)
inp$p3 <- c(1,1,0,2,1,0,2,0)
inp <- as.data.table(inp)
con <- list(hh1=35, hh2=65, p1=91, p2=65, p3=104) 
res <- ipu(inp=inp, hid="hhid", con=con, verbose=FALSE)
}
\keyword{method}

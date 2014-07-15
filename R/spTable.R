spTable <- function(inp, select) {
  if ( class(inp) != "synthPopObj") {
    stop("wrong input! Argument 'inp' must be of class 'synthPopObj'!\n")
  }
  dataS <- inp@sample@data
  dataP <- inp@pop@data
  weights <- inp@sample@weight
  
  # initializations
  # prepare weights
  if ( is.null(weights) ) {
    n <- nrow(dataS)
    weights <- rep.int(nrow(dataP)/n, n)
  } else {
    weights <- dataS[[weights]]
  }
  # prepare data seta
  if ( !all(select %in% colnames(dataS)) ) {
    stop("not all variables in argument 'select' available in the survey data available in argument 'inp'!\n")
  }
  if ( !all(select %in% colnames(dataP)) ) {
    stop("not all variables in argument 'select' available in the population data available in argument 'inp'!\n")
  }    
  dataS <- dataS[, select, with=FALSE]
  dataP <- dataP[, select, with=FALSE]

  # compute weighted table from sample (expected)
  tabS <- tableWt(dataS, weights)
  # compute table for population (realized)
  tabP <- table(dataP)
  # create object and return result
  res <- list(expected=tabS, realized=tabP)
  class(res) <- "spTable"
  res
}

# methods for class "spTable"
as.array.spTable <- function(x, ...) {
  values <- c(as.integer(x$expected), as.integer(x$realized))
  d <- c(dim(x$expected), 2)
  dnNew <- list(which=c("expected", "realized"))
  dn <- c(dimnames(x$expected), dnNew)
  array(values, dim=d, dimnames=dn)
}

as.table.spTable <- function(x, ...) {
  tab <- as.array(x)
  class(tab) <- "table"
  tab
}

plot.spTable <- function(x, ...) spMosaic(x, ...)

print.spTable <- function(x, ...) {
  # expected (from sample)
  cat("Expected:\n")
  print(x$expected, ...)
  # realized (population)
  cat("\nRealized:\n")
  print(x$realized, ...)
}

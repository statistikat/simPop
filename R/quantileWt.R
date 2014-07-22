quantileWt <- function(x, ...) UseMethod("quantileWt")

quantileWt.default <- function(x, weights=NULL, probs=seq(0, 1, 0.25), na.rm=TRUE, ...) {
  # initializations
  if ( !is.numeric(x) ) {
    stop("'x' must be a numeric vector!\n")
  }
  x <- unname(x)  # unlike 'quantile', this never returns a named vector
  if ( is.null(weights) ) {
    return(quantile(x, probs, na.rm=na.rm, names=FALSE, type=1))
  } else if ( !is.numeric(weights) ) {
    stop("'weights' must be a numeric vector!\n")
  } else if ( length(weights) != length(x) ) {
      stop("'weights' must have the same length as 'x'!\n")
  } else if ( !all(is.finite(weights)) ) {
    stop("missing or infinite weights!\n")
  }
  if ( !is.numeric(probs) || all(is.na(probs) ) || isTRUE(any(probs < 0 | probs > 1)) ) {
    stop("'probs' must be a numeric vector with values in [0,1]!\n")
  }
  if ( length(x) == 0 ) {
    return(rep.int(NA, length(probs)))
  }
  if ( !isTRUE(na.rm) && any(is.na(x)) ) {
    stop("missing values and NaN's not allowed if 'na.rm' is not TRUE")
  }
  # sort values and weights
  ord <- order(x, na.last=NA)
  x <- x[ord]
  weights <- weights[ord]
  # some preparations
  rw <- cumsum(weights)/sum(weights)
  # obtain quantiles
  select <- sapply(probs, function(p) min(which(rw >= p)))
  q <- x[select]
  invisible(q)
}

quantileWt.dataObj <- function(x, vars, probs=seq(0, 1, 0.25), na.rm=TRUE, ...) {
  dat <- x@data
  if ( is.null(dat) ) {
    return(NULL)
  } else {
    if ( length(vars) > 1 ) {
      stop("only one variable can be specified!\n")
    }
    ii <- match(vars, colnames(dat))
    if ( any(is.na(ii)) ) {
      stop("please provide valid variables that exist in the input object!\n")
    }
    tmpdat <- dat[[vars]]
    if ( !is.null(x@weight) ) {
      return(quantileWt.default(tmpdat, weights=dat[[x@weight]], probs=probs, na.rm=na.rm))
    } else {
      return(quantileWt.default(tmpdat, probs=probs, na.rm=na.rm))
    }
  }
}

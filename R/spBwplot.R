spBwplot <- function(inp, x, cond = NULL, horizontal = TRUE,
  coef = 1.5, zeros = TRUE, minRatio = NULL, do.out = FALSE, ...) {

  ## initializations
  if ( !class(inp) == "synthPopObj" ) {
    stop("input argument 'inp' must be of class 'synthPopObj'!\n")
  }

  weights.pop <- inp@pop@weight
  weights.samp <- inp@sample@weight
  dataS <- inp@sample@data
  dataP <- inp@pop@data

  if ( !is.character(x) || length(x) == 0 ) {
    stop("'x' must be a character vector of positive length!\n")
  }
  if ( !(all(x %in% colnames(dataP)) & (all(x %in% colnames(dataS)))) ) {
    stop("The variable names specified in argument 'x' must be available both in the population and the sample!\n")
  }

  if ( !is.null(cond) && !is.character(cond) ) {
    stop("'cond' must be a character vector or NULL!\n")
    if ( length(cond) != 1 ) {
      stop("argument 'cond' must have length 1!\n")
    }
  }
  if ( !(all(cond %in% colnames(dataP)) & (all(cond %in% colnames(dataS)))) ) {
    stop("The variable names specified in argument 'cond' must be available both in the population and the sample!")
  }

  horizontal <- isTRUE(horizontal)
  zeros <- isTRUE(zeros)
  do.out <- isTRUE(do.out)
  #nP <- if(inherits(dataP, "data.frame")) 2 else 1+length(dataP)
  nP <- 2
  pop <- "Population"
  lab <- c("Sample", pop)

  ## compute statistics for boxplots and construct objects for 'bwplot'
  # from sample
  tmp <- getBwplotStats(x, weights.samp, cond, dataS, coef=coef, zeros=zeros, do.out=do.out, name=lab[1])
  values <- tmp$values
  n <- t(tmp$n)
  nzero <- ifelse(zeros, t(tmp$nzero), NULL)
  out <- tmp$out
  # from population(s)
  tmp <- getBwplotStats(x, weights.pop, cond, dataP, coef=coef, zeros=zeros, do.out=do.out, name=lab[2])
  values <- rbind(values, tmp$values)
  n <- rbind(n, tmp$n)
  if ( zeros ) {
    nzero <- rbind(nzero, tmp$nzero)
  }
  out <- c(out, tmp$out)

  ## construct formula for 'bwplot'
  form <- ifelse(horizontal, ".name~.x", ".x~.name")  # basic formula
  if ( length(x) > 1 ) {
    cond <- c(".var", cond)
  }
  if ( !is.null(cond) ) {
    cond <- paste(cond, collapse = " + ")  # conditioning variabels
    form <- paste(form, cond, sep=" | ")  # add conditioning to formula
  }
  ## in case of semi-continuous variables define box widths
  if ( zeros ) {
    ratio <- n/(n+nzero)
    if ( !is.null(minRatio) ) {
      if ( !is.numeric(minRatio) || length(minRatio) != 1 || minRatio <= 0 || minRatio > 1 ) {
        stop("'minRatio' must be a single numeric value in [0,1]!\n")
      }
      ratio[ratio < minRatio] <- minRatio
    }
  } else {
    ratio <- NULL
  }
  ## define local version of 'bwplot'
  localBwplot <- function(form, values, xlab = NULL, ylab = NULL, ...,
  # these arguments are defined so that they aren't supplied twice:
  x, data, allow.multiple, outer, panel, groups) {
    bwplot(form, data=values, panel=panelSpBwplot, xlab=xlab, ylab=ylab, ...)
  }
  ## call 'bwplot'
  localBwplot(as.formula(form), values, horizontal=horizontal, coef=coef, zeros=zeros, ratio=ratio, do.out=FALSE, outliers=out, ...)
}

## panel function
panelSpBwplot <- function(x, y, coef=1.5, zeros = TRUE, ratio, outliers, subscripts, ...) {
  out <- outliers[subscripts]
  if ( zeros ) {
    i <- packet.number()
    localPanelBwplot <- function(..., ratio, box.ratio, box.width) {
      panel.bwplot(..., box.ratio=ratio)
    }
    localPanelBwplot(x[!out], y[!out], coef=coef, ratio=ratio[,i], ...)
  } else {
    panel.bwplot(x[!out], y[!out], coef=coef, ...)
  }
  panel.points(x[out], y[out], ...)
}


## internal utility functions

# get data.frame and all required statistics
getBwplotStats <- function(x, weights = NULL, cond = NULL, data, ..., name = "") {
  if ( is.null(cond) ) {
    x <- data[,x,with=FALSE]
    if ( length(weights) == 0 ) {
      w <- NULL
      } else {
        w <- data[[weights]]
    }
    prepBwplotStats(x, w, ..., name=name)
  } else {
    spl <- split(data, data[[cond]])
    tmp <- lapply(spl, function(z) {
      data <- z
      x <- data[,x,with=FALSE]
      if ( length(weights) == 0 ) {
        w <- NULL
      } else {
        w <- data[[weights]]
      }
      res <- prepBwplotStats(x, w, ..., name=name)
      res$values <- cbind(res$values, rep(data[[cond]][1], nrow(res$values)))
      colnames(res$values)[ncol(res$values)] <- cond
      res
    })
    values <- do.call(rbind, lapply(tmp, function(x) x$values))
    n <- as.vector(sapply(tmp, function(x) x$n))
    nzero <- as.vector(sapply(tmp, function(x) x$nzero))
    out <- do.call(c, lapply(tmp, function(x) x$out))
    list(values=values, n=n, nzero=nzero, out=out)
  }
}

# prepare one or more variables
prepBwplotStats <- function(x, w, ..., name = "") UseMethod("prepBwplotStats")

prepBwplotStats.data.frame <- function(x, w, ..., name = "") {
  tmp <- lapply(x, prepBwplotStats, w, ..., name=name)
  values <- mapply(function(x, v) cbind(x$values, .var=v), tmp, names(x), SIMPLIFY=FALSE, USE.NAMES=FALSE)
  values <- do.call(rbind, values)
  n <- sapply(tmp, function(x) x$n)
  nzero <- sapply(tmp, function(x) x$nzero)
  out <- do.call(c, lapply(tmp, function(x) x$out))
  list(values=values, n=n, nzero=nzero, out=out)
}

prepBwplotStats.default <- function(x, w, ..., name = "") {
  stats <- spBwplotStats(x, w, ...)
  x <- c(stats$stats, stats$out)
  n <- stats$n
  nzero <- stats$nzero
  nstats <- length(stats$stats)
  nout <- length(stats$out)
  out <- rep.int(c(FALSE, TRUE), c(nstats, nout))
  values <- data.frame(.x=x, .name=name)
  list(values=values, n=n, nzero=nzero, out=out)
}

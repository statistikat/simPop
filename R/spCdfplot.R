spCdfplot <- function(inp, x, cond = NULL, approx = c(FALSE, TRUE), 
                      n = 10000, bounds = TRUE, ...) {
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
  
  # check 'approx'
  if(!is.logical(approx) || length(approx) == 0) approx <- formals()$approx
  else approx <- sapply(rep(approx, length.out=2), isTRUE)
  # check 'n'
  if(any(approx) && (!is.numeric(n) || length(n) == 0)) {
    stop("'n' is not numeric or does not have positive length")
  } else n <- ifelse(approx, n[1], NA)
  # check 'bounds'
  bounds <- isTRUE(bounds)
  # define labels for grouping variable
  lab <- c("Sample", "Population")
  
  ## construct objects for 'xyplot'
  # from sample
  tmp <- getCdf(x, weights.samp, cond, dataS, approx=approx[1], n=n[1], name=lab[1])
  values <- tmp$values
  app <- t(tmp$approx)
  # from population
  tmp <- getCdf(x, weights.pop, cond, dataP, approx=approx[2], n=n[2], name=lab[2])
  values <- rbind(values, tmp$values)
  app <- rbind(app, tmp$approx)
  ## construct formula for 'xyplot'
  form <- ".y~.x"  # basic formula
  if(length(x) > 1) cond <- c(".var", cond)
  if(!is.null(cond)) {
    cond <- paste(cond, collapse = " + ")  # conditioning variabels
    form <- paste(form, cond, sep=" | ")  # add conditioning to formula
  }
  ## define local version of 'xyplot'
  localXyplot <- function(form, values, xlab = NULL, ylab = NULL, 
                          auto.key = TRUE, ..., 
                          # these arguments are defined so that they aren't supplied twice:
                          x, data, allow.multiple, outer, panel, prepanel, groups) {
    # prepare legend
    if(isTRUE(auto.key)) auto.key <- list(points=FALSE, lines=TRUE)
    else if(is.list(auto.key)) {
      if(is.null(auto.key$points)) auto.key$points <- FALSE
      if(is.null(auto.key$lines)) auto.key$lines <- TRUE
    }
    command <- paste("xyplot(form, data=values, groups=.name,", 
                     "panel=panelSpCdfplot, prepanel=prepanelSpCdfplot,", 
                     "xlab=xlab, ylab=ylab, auto.key=auto.key, ...)")
    eval(parse(text=command))
  }
  ## call 'xyplot'
  localXyplot(as.formula(form), values, approx=app, bounds=bounds, ...)
}


## panel function
panelSpCdfplot <- function(x, y, approx, bounds = TRUE, ...) {
  if(isTRUE(bounds)) {
    panel.refline(h=0, ...)
    panel.refline(h=1, ...)
  }
  localPanelXyplot <- function(..., approx, type, distribute.type) {
    i <- packet.number()
    type <- ifelse(approx[,i], "l", "s")
    panel.xyplot(..., type=type, distribute.type=TRUE)
  }
  localPanelXyplot(x, y, approx=approx, ...)
}

## prepanel function
prepanelSpCdfplot <- function(x, y, ...) list(ylim=c(0,1))


## internal utility functions

# get data.frame and logical indicating approximation
getCdf <- function(x, weights = NULL, cond = NULL, data, ..., name = "") {
  if ( is.null(cond) ) {
    x <- data[, x, with=FALSE]
    if ( length(weights) == 0 ) {
      w <- NULL
    } else {
      w <- data[[weights]]
    }
    prepCdf(x, w, ..., name=name)
  } else {
    spl <- split(data, data[[cond]])
    tmp <- lapply(spl, function(z) {
      data <- z
      x <- data[, x, with=FALSE]
      if ( length(weights) == 0 ) {
        w <- NULL
      } else {
        w <- data[[weights]]
      }
      res <- prepCdf(x, w, ..., name=name)
      res$values <- cbind(res$values, rep(data[[cond]][1], nrow(res$values)))
      colnames(res$values)[ncol(res$values)] <- cond
      res
    })
    values <- do.call(rbind, lapply(tmp, function(x) x$values))
    approx <- as.vector(sapply(tmp, function(x) x$approx))
    list(values=values, approx=approx)
  }
}

# prepare one or more variables
prepCdf <- function(x, w, ..., name = "") UseMethod("prepCdf")

prepCdf.data.frame <- function(x, w, ..., name = "") {
  tmp <- lapply(x, prepCdf, w, ..., name=name)
  values <- mapply(function(x, v) cbind(x$values, .var=v), tmp, names(x), SIMPLIFY=FALSE, USE.NAMES=FALSE)
  values <- do.call(rbind, values)
  approx <- sapply(tmp, function(x) x$approx)
  list(values=values, approx=approx)
}

prepCdf.default <- function(x, w, ..., name = "") {
  tmp <- spCdf(x, w, ...)
  values <- data.frame(.x=c(tmp$x[1], tmp$x), .y=c(0, tmp$y), .name=name)
  list(values=values, approx=tmp$approx)
}

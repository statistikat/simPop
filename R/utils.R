# calculate parameters for parallel processing
# returns a list with following elements
# - have_win: TRUE if windows system, FALSE else
# - nr_cores: number of cpus to use
# - parallel: TRUE if parallel computing should be applied
parallelParameters <- function(nr_cpus=NULL, nr_strata) {
  control <- list()
  control$have_win <- Sys.info()["sysname"] == "Windows"

  cpus_available <- detectCores()
  if ( !is.null(nr_cpus) && nr_cpus > cpus_available) {
    stop("there are only",cpus_available,"cpus available in your system!\n")
  }
  if ( !is.null(nr_cpus) && nr_cpus < 1) {
    stop("we must use at least one cpu!\n")
  }

  if ( !is.null(nr_cpus) ) {
    control$nr_cores <- min(nr_cpus, nr_strata)
  } else {
    if ( cpus_available > 1 & nr_strata > 1 ) {
      control$nr_cores <- min(cpus_available-1, nr_strata) # keep at least one core available
    } else {
      control$nr_cores <- 1
    }
  }
  control$parallel <- ifelse(control$nr_cores>1, TRUE, FALSE)
  control
}

## check break points for categorization
checkBreaks <- function(x) {
  if(!is.numeric(x) || length(x) < 2) {
    stop("'breaks' must be a numeric vector with length >= 2")
  }
  invisible()
}

## check whether selected columns of a data frame are factors
# TODO: add argument 'num': a logical indicating whether numeric variables other than integers should be converted to factors
checkFactor <- function(x, select) {
  if ( "data.table" %in% class(x) ) {
    convert <- select[sapply(x[, select, with=F], function(x) !is.factor(x))]
    if ( length(convert) > 0 ) {
      x[,(convert):=lapply(.SD, as.factor),.SDcols=convert]
    }
  } else {
    convert <- select[sapply(x[, select], function(x) !is.factor(x))]
    n <- length(convert)
    if ( n == 1 ) {
      x[[convert]] <- as.factor(x[[convert]])
    } else {
      x[, convert] <- as.data.frame(lapply(x[, convert], as.factor))
    }
  }
  invisible(x)
}

# drop unused factor levels for factors in a data.frame
dropLevels <- function(x, select = names(x)) {
  check <- select[sapply(x[, select], function(x) is.factor(x))]
  n <- length(check)
  if(n == 1) x[, check] <- x[, check][, drop=TRUE]
  else if(n > 1) {
    x[, check] <- as.data.frame(lapply(x[, check], "[", drop=TRUE))
  }
  x
}

# create factor that includes NA
factorNA <- function(x, always = FALSE) {
  always <- isTRUE(always)
  if(is.factor(x)) {
    l <- levels(x)
    if(NA %in% l || !(always || any(is.na(x)))) x
    else {
      l <- c(l, NA)
      factor(x, levels=c(levels(x), NA), exclude=c())
    }
  } else {
    if(always) {
      factor(c(NA, x), exclude=c())[-1]  # little trick
    } else factor(x, exclude=c())
  }
}


## get which observations contain NAs (and need to be excluded)
getExclude <- function(x, ...) UseMethod("getExclude")
getExclude.default <- function(x, ...) which(is.na(x))
getExclude.data.frame <- function(x, ...) {
  unique(which(is.na(x), arr.ind=TRUE)[, 1])
}
getExclude.data.table <- function(x, ...) {
  unique(which(is.na(x), arr.ind=TRUE)[, 1])
}


### exclude observations
#excludeData <- function(x, exclude = NULL, ...) UseMethod("excludeData")
#excludeData.default <- function(x, exclude = NULL, ...) {
#	if(length(exclude) == 0) x else x[-exclude]
#}
#excludeData.data.frame <- function(x, exclude = NULL, ...) {
#	if(length(exclude) == 0) x else x[-exclude, , drop=FALSE]
#}

## adjust probabilities estimated with multinomial model to account for
## structural zeros
adjustProbs <- function(probs, grid, pNames, limit = NULL, censor = NULL) {
  set0 <- matrix(FALSE, nrow(probs), ncol(probs),
                 dimnames=dimnames(probs))
  target <- colnames(probs)
  # account for structural zeros via argument 'limit'
  if(is.list(limit)) {
    limit <- limit[sapply(limit, inherits, "list")]
    limit <- limit[names(limit) %in% names(grid)]
    # loop over predictors for which to censor probabilities
    for(i in seq_along(limit)) {
      predNameI <- names(limit)[i]
      predI <- grid[, predNameI]
      limitI <- limit[[i]]
      # loop over supplied outcomes of current predictor to
      # find probabilities to be set to zero
      for(j in seq_along(limitI)) {
        cat <- names(limitI)[j]
        ok <- intersect(target, limitI[[j]])
        if(length(ok) > 0) {
          set0[predI == cat, setdiff(target, ok)] <- TRUE
        }
      }
    }
  }
  # account for structural zeros via argument 'censor'
  if(is.list(censor)) {
    censor <- censor[sapply(censor, inherits, c("list", "data.frame"))]
    censor <- censor[names(censor) %in% target]
    if(length(censor)) {
      for(i in seq_along(censor)) {
        censorNameI <- names(censor)[i]
        censorI <- censor[[i]]
        if(is.list(censorI)) {
          # loop over supplied predictors to find
          # probabilities to be set to zero
          for(j in seq_along(censorI)) {
            predNameJ <- names(censorI)[j]
            predJ <- grid[, predNameJ]
            set0[predJ %in% censorI[[j]], censorNameI] <- TRUE
          }
        } else {
          # set probabilities to zero for supplied
          # combinations of predictors
          cNames <- apply(censorI, 1, paste, collapse=".")
          set0[pNames %in% cNames, censorNameI] <- TRUE
        }
      }
    }
  }
  # set indicated probabilities to zero
  probs[set0] <- 0
  # set the non-censored probabilities to non-zero if all
  # probabilities of a row are zero
  adjust <- which(apply(probs == 0, 1, all))
  for(i in adjust) {
    ok <- which(!set0[i,])
    probs[i, ok] <- 1/length(ok)
  }
  probs
}

## get breakpoints for categorizing continuous or semi-continuous variables
getBreaks <- function(x, weights = NULL, zeros = TRUE, lower = NULL,
                      upper = NULL, equidist = TRUE, probs = NULL) {
  # initializations
  if(!is.numeric(x)) stop("'x' must be a numeric vector")
  if(!is.null(weights)) {
    if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
    else if(length(weights) != length(x)) {
      stop("'weights' must have the same length as 'x'")
    }
  }
  zeros <- isTRUE(zeros)
  if(!is.null(probs)) {
    if(!is.numeric(probs) || all(is.na(probs)) ||
       isTRUE(any(probs < 0 | probs > 1))) {
      stop("'probs' must be a numeric vector with values in [0,1]")
    }
  }
  if(zeros) {
    pos <- which(x > 0)
    if(length(pos)) {
      if(is.null(probs)) {
        if(isTRUE(equidist)) probs <- seq(0.1, 1, by=0.1)
        else probs <- c(0.01, 0.05, 0.1, seq(0.2, 0.8, by=0.2), 0.9, 0.95, 0.99, 1)
      } else probs <- c(probs, 1)
      qpos <- quantileWt(x[pos], weights[pos], probs)
    } else qpos <- NULL
    neg <- which(x < 0)
    if(length(neg)) {
      pneg <- seq(0, 0.9, by=0.1)
      qneg <- quantileWt(x[neg], weights[neg], pneg)
    } else qneg <- NULL
    breaks <- c(qneg, 0, qpos)
  } else {
    if(is.null(probs)) {
      if(isTRUE(equidist)) probs <- seq(0, 1, by=0.1)
      else probs <- c(0, 0.01, 0.05, 0.1, seq(0.2, 0.8, by=0.2), 0.9, 0.95, 0.99, 1)
    } else probs <- c(0, probs, 1)
    breaks <- quantileWt(x, weights, probs)
  }
  breaks <- unique(breaks)  # remove duplicated values
  if(!is.null(lower)) {
    if(!is.numeric(lower) || length(lower) == 0) {
      stop("'lower' is not numeric or has length 0")
    } else if(length(lower) > 1) lower <- lower[1]
    if(isTRUE(lower > breaks[1])) {
      warning("'lower' is larger than the smallest ",
              "breakpoint and therefore disregarded")
    } else breaks[1] <- lower
  }
  if(!is.null(upper)) {
    if(!is.numeric(upper) || length(upper) == 0) {
      stop("'upper' is not numeric or has length 0")
    } else if(length(upper) > 1) upper <- upper[1]
    nb <- length(breaks)
    if(isTRUE(upper < breaks[nb])) {
      warning("'upper' is smaller than the largest ",
              "breakpoint and therefore disregarded")
    } else breaks[nb] <- upper
  }
  breaks
}


## categorize continuous or semi-continuous variables
getCat <- function(x, breaks, zeros = TRUE, right = FALSE) {
  # initializations
  if(!is.numeric(x)) stop("'x' must be a numeric vector")
  checkBreaks(breaks)
  if(isTRUE(zeros)) {
    pos <- which(x > 0)
    neg <- which(x < 0)
    # positive values (also works if none exist)
    bpos <- c(0, breaks[breaks > 0])
    if(length(bpos) == 1) lpos <- NULL
    else {
      cpos <- cut(x[pos], bpos)
      lpos <- levels(cpos)
    }
    # negative values (also works if none exist)
    bneg <- c(breaks[breaks < 0], 0)
    if(length(bneg) == 1) lneg <- NULL
    else {
      cneg <- cut(x[neg], bneg, right=FALSE)
      lneg <- levels(cneg)
    }
    # put it all together
    categories <- factor(ifelse(is.na(x), NA, 0), levels=c(lneg, 0, lpos))
    if(length(pos) > 0 && length(bpos) > 1) categories[pos] <- cpos
    if(length(neg) > 0 && length(bneg) > 1) categories[neg] <- cneg
    # return vector
    categories
  } else cut(x, breaks, include.lowest=TRUE, right=right)
}


## get name of categorized variable
getCatName <- function(name) paste(name, "Cat", sep="")


## get name of categorical variable for household head
getHeadName <- function(name) paste(name, "Head", sep="")


## truncated Pareto distribution
truncPareto <- function(n, loc, scale, shape, lower, upper) {
  # initializations
  x <- numeric(n)
  left <- n
  # fit generalized pareto distribution with lower and upper bound
  while(left > 0) {
    xTmp <- rgpd(left, loc=loc, scale=scale, shape=shape)
    ind  <- which(xTmp > lower & xTmp <= upper)
    lind <- length(ind)
    if(lind > 0) {
      x[(n - left + 1):(n - left + lind)]  <- xTmp[ind]
      left <- left - lind
    }
  }
  return(x)
}


## logit regression (designed for internal use, hence no error handling)
#logitreg <- function(x, y, weights = rep(1, length(y)),
#    intercept = TRUE, start = rep(0, p), ...) {
#    # function to be minimized (log-likelihood)
#    fmin <- function(beta, X, y, w) {
#        p <- plogis(X %*% beta)
#        -sum(2 * w * ifelse(y, log(p), log(1-p)))
#    }
#    # gradient
#    gmin <- function(beta, X, y, w) {
#        eta <- as.numeric(X %*% beta)
#        p <- plogis(eta)
#        -2 * (w * dlogis(eta) * ifelse(y, 1/p, -1/(1-p))) %*% X
#    }
#    # some preparations
#    if(is.null(dim(x))) dim(x) <- c(length(x), 1)
#    dn <- dimnames(x)[[2]]
#    if(!length(dn)) dn <- paste("Var", 1:ncol(x), sep="")
#    p <- ncol(x) + intercept
#    if(intercept) {
#        x <- cbind(1,x)
#        dn <- c("(Intercept)", dn)
#    }
#    if(is.factor(y)) y <- (unclass(y) != 1)
#    # optimize and return result
#    fit <- optim(start, fmin, gmin, X = x, y = y, w = weights, method = "BFGS", ...)
#    names(fit$par) <- dn
#    return(fit)
#}
# intercept is assumed to be included in 'x'
logitreg <- function(x, y, weights = rep(1, length(y)), start = rep(0, p), ...) {
  # function to be minimized (log-likelihood)
  fmin <- function(beta, X, y, w) {
    p <- plogis(X %*% beta)
    -sum(2 * w * ifelse(y, log(p), log(1-p)))
  }
  # gradient
  gmin <- function(beta, X, y, w) {
    eta <- as.numeric(X %*% beta)
    p <- plogis(eta)
    -2 * (w * dlogis(eta) * ifelse(y, 1/p, -1/(1-p))) %*% X
  }
  # some preparations
  if(is.null(dim(x))) dim(x) <- c(length(x), 1)
  dn <- dimnames(x)[[2]]
  if(!length(dn)) dn <- paste("Var", 1:ncol(x), sep="")
  p <- ncol(x)
  if(is.factor(y)) y <- (unclass(y) != 1)
  # optimize and return result
  fit <- optim(start, fmin, gmin, X = x, y = y, w = weights, method = "BFGS", ...)
  names(fit$par) <- dn
  return(fit)
}


## search for donor cells for splitting continuous variables
# i ........... indices corresponding to a zero of a table
# indDonors ... matrix containing indices of all non-zero elements
# donors ...... vector indices of the non-zero elements

# sequential minimum distance
seqMinDist <- function(i, indDonors, donors) {
  # initializations
  donor <- NA
  k <- 1
  # loop over dimensions to find donor cells with matching category
  while(is.na(donor) && k <= length(i)) {
    sel <- which(indDonors[, k] == i[k])  # donor cells in same category
    if(length(sel)) {
      tmpIndDonors <- indDonors[sel, , drop=FALSE]
      tmpDonors <- donors[sel]
      donor <- minDist(i[-k], tmpIndDonors[, -k, drop=FALSE], tmpDonors)
    }
    k <- k + 1
  }
  # if not successful, use donor with minimum overall distance
  if(is.na(donor)) donor <- minDist(i, indDonors, donors)
  donor
}

# minimum distance
minDist <- function(i, indDonors, donors) {
  # calculate distance from element defined by i for each element of m
  d <- colSums(abs(t(indDonors) - i))
  # return index of selected donor (minimum distance)
  donors[which.min(d)]
}


## weighted mean
meanWt <- function(x, ...) UseMethod("meanWt")

meanWt.default <- function(x, weights, na.rm=TRUE, ...) {
  na.rm <- isTRUE(na.rm)
  if(missing(weights)) mean(x, na.rm=na.rm)
  else weighted.mean(x, w=weights, na.rm=na.rm)
}

meanWt.dataObj <- function(x, vars, na.rm=TRUE, ...) {
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
      return(meanWt.default(tmpdat, weights=dat[[x@weight]], na.rm=na.rm))
    } else {
      return(meanWt.default(tmpdat, na.rm=na.rm))
    }
  }
}


## weighted variance
varWt <- function(x, ...) UseMethod("varWt")

varWt.default <- function(x, weights, na.rm=TRUE, ...) {
  na.rm <- isTRUE(na.rm)
  if(missing(weights)) var(x, na.rm=na.rm)
  else {
    x <- as.numeric(x)
    weights <- as.numeric(weights)
    if(length(weights) != length(x)) {
      stop("'weights' must have the same length as 'x'")
    }
    if(na.rm) {
      select <- !is.na(x)
      x <- x[select]
      weights <- weights[select]
    }
    if(length(x) <= 1 || sum(weights > 0) <= 1) NA
    else sum((x - meanWt(x, weights))^2 * weights) / (sum(weights) - 1)
  }
}

varWt.dataObj <- function(x, vars, na.rm=TRUE, ...) {
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
      return(varWt.default(tmpdat, weights=dat[[x@weight]], na.rm=na.rm))
    } else {
      return(varWt.default(tmpdat, na.rm=na.rm))
    }
  }
}

## weighted covariance matrix

# generic function
covWt <- function(x, ...) UseMethod("covWt")

# default method
covWt.default <- function(x, y, weights, ...) {
  if(missing(y)) y <- x
  else if(length(x) != length(y)) {
    stop("'x' and 'y' must have the same length")
  }
  if(missing(weights)) cov(x, y, use = "complete.obs")
  else {
    if(length(weights) != length(x)) {
      stop("'weights' must have the same length as 'x' and 'y'")
    }
    select  <- !is.na(x) & !is.na(y)
    x <- x[select]
    y <- y[select]
    weights <- weights[select]
    sum((x-meanWt(x, weights)) * (y-meanWt(y, weights)) * weights) /
    (sum(weights)-1)
  }
}

# method for matrices
covWt.matrix <- function(x, weights, ...) {
  if(missing(weights)) cov(x, use = "pairwise.complete.obs")
  else {
    if(length(weights) != nrow(x)) {
      stop("length of 'weights' must equal the number of rows in 'x'")
    }
    center <- apply(x, 2, meanWt, weights=weights)
    x <- sweep(x, 2, center, check.margin = FALSE)
    crossprodWt(x, weights)
  }
}

# method for data.frames
covWt.data.frame <- function(x, weights, ...) covWt(as.matrix(x), weights)

covWt.dataObj <- function(x, vars, ...) {
  dat <- x@data
  if ( is.null(dat) ) {
    return(NULL)
  } else {
    ii <- match(vars, colnames(dat))
    if ( any(is.na(ii)) ) {
      stop("please provide valid variables that exist in the input object!\n")
    }
    tmpdat <- dat[,vars,with=F]
    if ( !is.null(x@weight) ) {
      return(covWt.matrix(as.matrix(tmpdat), weights=dat[[x@weight]]))
    } else {
      return(covWt.matrix(as.matrix(tmpdat)))
    }
  }
}

## weighted correlation matrix

# generic function
corWt <- function(x, ...) UseMethod("corWt")

# default method
corWt.default <- function(x, y, weights, ...) {
  if(missing(y)) y <- x
  else if(length(x) != length(y)) {
    stop("'x' and 'y' must have the same length")
  }
  if(missing(weights)) cor(x, y, use = "complete.obs")
  else covWt(x, y, weights) / sqrt(varWt(x, weights) * varWt(y, weights))
}

# method for matrices
corWt.matrix <- function(x, weights, ...) {
  if(missing(weights)) cor(x, use = "pairwise.complete.obs")
  else {
    if(length(weights) != nrow(x)) {
      stop("length of 'weights' must equal the number of rows in 'x'")
    }
    cen <- apply(x, 2, meanWt, weights=weights)
    sc <- sqrt(apply(x, 2, varWt, weights=weights))
    x <- scale(x, center=cen, scale=sc)
    crossprodWt(x, weights)
  }
}

# method for data.frames
corWt.data.frame <- function(x, weights, ...) corWt(as.matrix(x), weights)

# method for objects of class "dataObj"
corWt.dataObj <- function(x, vars, ...) {
  dat <- x@data
  if ( is.null(dat) ) {
    return(NULL)
  } else {
    ii <- match(vars, colnames(dat))
    if ( any(is.na(ii)) ) {
      stop("please provide valid variables that exist in the input object!\n")
    }
    tmpdat <- dat[,vars,with=F]
    if ( !is.null(x@weight) ) {
      return(corWt.matrix(as.matrix(tmpdat), weights=dat[[x@weight]]))
    } else {
      return(corWt.matrix(as.matrix(tmpdat)))
    }
  }
}

## weighted cross product
# designed for internal use, hence no error handling
# TODO: more efficient solution
crossprodWt <- function(x, weights) {
  ci <- 1:ncol(x)
  sapply(ci, function(j) sapply(ci, function(i) {
    select <- !is.na(x[, i]) & !is.na(x[, j])
    xi <- x[select, i]
    xj <- x[select, j]
    w <- weights[select]
    sum(xi*xj*w) / (sum(w)-1)
  }))
}

manageSimPopObj <- function(x, var, sample=FALSE, set=FALSE, values=NULL) {
  if ( class(x) != "simPopObj" ) {
    stop("wrong input of argument 'x' (needs to be of class 'simPopObj')!\n")
  }
  if ( length(var) != 1 ) {
    stop("only one variable can be used at a time!\n")
  }
  if ( set==FALSE ) {
    if ( sample ) {
      return(invisible(samp(x, var=var)))
    } else {
      return(invisible(pop(x, var=var)))
    }
  }
  if ( set == TRUE ) {
    if ( is.null(values) ) {
      stop("you need to provide values!\n")
    }
    if ( sample ) {
      samp(x, var=var) <- values
    } else {
      pop(x, var=var) <- values
    }
    return(invisible(x))
  }
}

setGeneric("samp", function(obj, var) {
  standardGeneric("samp")
})
setMethod("samp", "simPopObj", function(obj, var) {
  if ( is.numeric(var) ) {
    if ( !all(var <= ncol(obj@sample@data)) ) {
      stop("check input 'var'!\n")
    }
  }
  if ( is.character(var) ) {
    cn <- colnames(obj@sample@data)
    var <- var[var%in%cn]
    if ( length(var) == 0 ) {
      return(NULL)
    }
  }
  if ( is.data.table(obj@pop@data)) {
    obj@sample@data[,var,with=F]
  } else {
    obj@sample@data[,var,drop=F]
  }
})
setGeneric("samp<-", function(obj, var, value) {
  standardGeneric("samp<-")
})
setReplaceMethod("samp", "simPopObj", function(obj, var, value) {
  if ( length(var) != 1) {
    stop("we can only set one variable!\n")
  }
  obj@sample@data[[var]] <- value
  validObject(obj)
  obj
})

setGeneric("pop", function(obj, var) {
  standardGeneric("pop")
})
setMethod("pop", "simPopObj", function(obj, var) {
  if ( is.numeric(var) ) {
    if ( !all(var <= ncol(obj@pop@data)) ) {
      stop("check input 'var'!\n")
    }
  }
  if ( is.character(var) ) {
    cn <- colnames(obj@pop@data)
    var <- var[var%in%cn]
    if ( length(var) == 0 ) {
      return(NULL)
    }
  }
  if ( is.data.table(obj@pop@data)) {
    obj@pop@data[,var,with=F]
  } else {
    obj@pop@data[,var,drop=F]
  }
})
setGeneric("pop<-", function(obj, var, value) {
  standardGeneric("pop<-")
})
setReplaceMethod("pop", "simPopObj", function(obj, var, value) {
  if ( length(var) != 1) {
    stop("we can only set one variable!\n")
  }
  obj@pop@data[[var]] <- value
  validObject(obj)
  obj
})

################################################################
### utility functions for the sga project
getAge <- function(birth, year, data = NULL) {
  # this function is only applicable if the income 
  # reference period is the previous calender year
  if(!is.null(data)) {
    if(missing(birth)) birth <- "rb080"
    birth <- data[, birth]
    if(missing(year)) year <- "rb010"
    year <- data[, year]
  }
  if(is.factor(year)) year <- as.numeric(as.character(year))
  if(is.factor(birth)) year <- as.numeric(as.character(birth))
  year - 1 - birth
}
getGender <- function(gender, labels = c("male","female"), data = NULL) {
  if(!is.null(data)) {
    if(missing(gender)) gender <- "rb090"
    gender <- c(data[, gender])
  }
  if(is.list(gender)) gender <- gender[[1]]
  factor(c(gender), labels=labels)
}
#getHsize <- function(data) 
#{
#  tab <- table(data$rb040)
#  hsize <- rep(tab, tab)
#}
getEcoStat <- function(ecoStat , data , levels) {  ## variable pl030 (economic status)
  if(missing(ecoStat)) ecoStat <- "pl031"
  ecoStat <- factor(data[, ecoStat])
  levels(ecoStat) <- levels
  return(ecoStat)
}
getCitizenship <- function(citizenship, data, owncountry, EU, other) {
  if(missing(citizenship)) citizenship <- "pb220a"
  citizenship <- data[, citizenship]
  indNA <- which(levels(citizenship) == "")
  indOC <- which(levels(citizenship) == owncountry)
  indEU <- which(levels(citizenship) %in% EU)
  indOther <- which(levels(citizenship) %in% other)
  levels <- character(nlevels(citizenship))
  levels[indNA] <- NA
  levels[indOC] <- owncountry
  levels[indEU] <- "EU"
  levels[indOther] <- "Other"
  levels(citizenship) <- levels
  return(citizenship)
}
getHsize <- function(data,hhid) 
{ 
  if(missing(hhid)) hhid <- "db030"
  tab <- table(data[,hhid]) #table(data$rb040)
  hsize <- rep(tab, tab)
  hsize <- as.numeric((hsize))
  return(hsize)
}
restructureHHid <- function(data){
  tab <- table(data$db030)
  hsize <- rep(tab, tab)
  db030 <- as.numeric(names(hsize))
  return(db030)
}
#Function factorNA from package simPop: includes NAs as an extra level in the factor
factorNA <- function(x, always = FALSE) {
  always <- isTRUE(always)
  if(is.factor(x)) {
    l <- levels(x)
    if(NA %in% l || !(always || any(is.na(x)))) x
    else {
      l <- c(l, NA)
      factor(x, levels=c(levels(x), NA), exclude=c())
    }
  } else {
    if(always) {
      factor(c(NA, x), exclude=c())[-1] # little trick
    } else factor(x, exclude=c())
  }
}

# Function uni.distribution: random draws from the weighted univariate distribution of
# the original data (maybe better from the SUF, but then the SUF always has to be used as well)
univariate.dis <- function(puf,data,additional,w){
  if (sum(is.na(data[,additional]))>0 & sum(is.na(data[,additional])) != dim(data)[1]) {
    var <- factorNA(data[,additional],always=TRUE)
  } else if (sum(is.na(data[,additional])) == dim(data)[1]) {
    var <- factor(c(NA, data[,additional]), exclude=c())[-1]
  } else {
    var <- as.factor(data[,additional])
  }
  tab <- wtd.table(var,weights=data[,w],type="table")
  p <- tab/sum(data[,w])
  puf[,additional] <- sample(x=levels(var)[levels(var) %in% names(tab)],size=dim(puf)[1],prob=p,replace=T)
  return(puf)
}

# Function con.distribution: random draws from the weighted conditional distribution
# (conditioned on a factor variable)
conditional.dis <- function(puf,data,additional,conditional,w){
  if (sum(is.na(data[,additional]))>0 & sum(is.na(data[,additional])) != dim(data)[1]) {
    var <- factorNA(data[,additional],always=TRUE)
  } else if (sum(is.na(data[,additional])) == dim(data)[1]) {
    var <- factor(c(NA, data[,additional]), exclude=c())[-1]
  } else {
    var <- as.factor(data[,additional])
  }
  puf[,additional] <- NA
  for (i in 1:length(levels(puf[,conditional]))) {
    tab <- wtd.table(var[data[,conditional]==levels(data[,conditional])[i]],weights=data[data[,conditional]==levels(data[,conditional])[i],w],type="table")
    p <- tab/sum(tab)
    puf[which(puf[,conditional]==levels(puf[,conditional])[i]),additional] <- sample(x=levels(var)[levels(var) %in% names(tab)],size=dim(puf[which(puf[,conditional]==levels(data[,conditional])[i]),])[1],prob=p,replace=T)
  }
  return(puf)
}
utility <- function(x, y, type="all"){
  if(type=="all" | type=="measure2"){
    measure2 <- ncol(x) / ncol(y)
  }
  if(type=="all" | type=="measure3"){
    measure3 <- nrow(x) / nrow(y)
  }
  if(type=="all" | type=="measure4"){
    puf <- sum(is.na(x))
    suf <- sum(is.na(y))
    if(suf > 0 & puf > 0){ 
      measure4 <- sum(is.na(x)) / sum(is.na(y)) - 1
    } else if(suf == 0 & puf == 0) {
      measure4 <- 0  
    } else if(suf > 0 & puf == 0){
      measure4 <- 1
    } else if(suf == 0 & puf > 0){
      measure4 <- min(c(puf / nrow(x), 1))
    }
  }
  measures <- list("measure2"=measure2,
                   "measure3"=measure3,
                   "measure4"=measure4)
  return(measures)
}
utilityModal <- function(x, y, variable){
  measure5 <- length(table(x[, variable])) / length(table(y[, variable]))
  return(measure5)
}
utilityIndicator <- function(indicatorPUF, indicatorSUF){
  measure6 <- abs(indicatorPUF - indicatorSUF) / indicatorSUF 
  return(measure6)
}
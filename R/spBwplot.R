# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

spBwplot <- function(x, ...) UseMethod("spBwplot")

spBwplot.default <- function(x, weights = NULL, cond = NULL, dataS, 
        dataP = NULL, horizontal = TRUE, coef = 1.5, zeros = TRUE, 
        minRatio = NULL, do.out = FALSE, ...) {
    ## initializations
    if(!is.character(x) || length(x) == 0) {
        stop("'x' must be a character vector of positive length")
    }
    if(!is.null(weights) && (!is.character(weights) || length(weights) > 1)) {
        stop("'weights' must be a single character string or NULL")
    }
    if(!is.null(cond) && !is.character(cond)) {
        stop("'cond' must be a character vector or NULL")
    }
    if(!inherits(dataS, "data.frame")) stop("'dataS' must be a data.frame")
    # check if 'dataP' is valid
    if(!is.null(dataP) && !inherits(dataP, "data.frame")) {
        if(ok <- inherits(dataP, "list")) {
            if(length(dataP)) ok <- all(sapply(dataP, inherits, "data.frame"))
            else dataP <- NULL
        }
        if(!ok) stop("'data' must be a data.frame or a list of data.frames")
    }
    horizontal <- isTRUE(horizontal)
    zeros <- isTRUE(zeros)
    do.out <- isTRUE(do.out)
    nP <- if(inherits(dataP, "data.frame")) 2 else 1+length(dataP)
    if(is.null(dataP)) lab <- ""
    else {
        pop <- "Population"
        if(!inherits(dataP, "data.frame")) {
            nam <- names(dataP)
            if(is.null(nam)) pop <- paste(pop, 1:(nP-1))
            else {
                replace <- which(nchar(nam) == 0)
                nam[replace] <- paste(pop, replace)
                pop <- nam
            }
        }
        lab <- c("Sample", pop)
    }
    ## compute statistics for boxplots and construct objects for 'bwplot'
    # from sample
    tmp <- getBwplotStats(x, weights, cond, dataS, 
        coef=coef, zeros=zeros, do.out=do.out, name=lab[1])
    values <- tmp$values
    n <- t(tmp$n)
    nzero <- if(zeros) t(tmp$nzero) else NULL
    out <- tmp$out
    # from population(s)
    if(!is.null(dataP)) {
        if(inherits(dataP, "data.frame")) {
            tmp <- getBwplotStats(x, NULL, cond, dataP, 
                coef=coef, zeros=zeros, do.out=do.out, name=lab[2])
            values <- rbind(values, tmp$values)
            n <- rbind(n, tmp$n)
            if(zeros) nzero <- rbind(nzero, tmp$nzero)
            out <- c(out, tmp$out)
        } else {
            tmp <- mapply(function(dP, l) {
                    getBwplotStats(x, NULL, cond, dP, coef=coef, 
                        zeros=zeros, do.out=do.out, name=l)
                }, dataP, lab[-1], SIMPLIFY=FALSE, USE.NAMES=FALSE)
            values <- rbind(values, 
                do.call(rbind, lapply(tmp, function(x) x$values)))
            n <- rbind(n, do.call(rbind, lapply(tmp, function(x) x$n)))
            if(zeros) {
                nzero <- rbind(nzero, 
                    do.call(rbind, lapply(tmp, function(x) x$nzero)))
            }
            out <- c(out, do.call(c, lapply(tmp, function(x) x$out)))
        }
    }
    ## construct formula for 'bwplot'
    form <- if(horizontal) ".name~.x" else ".x~.name"  # basic formula
    if(length(x) > 1) cond <- c(".var", cond)
    if(!is.null(cond)) {
        cond <- paste(cond, collapse = " + ")  # conditioning variabels
        form <- paste(form, cond, sep=" | ")  # add conditioning to formula
    }
    ## in case of semi-continuous variables define box widths 
    if(zeros) {
        ratio <- n/(n+nzero)
        if(!is.null(minRatio)) {
            if(!is.numeric(minRatio) || length(minRatio) != 1 || 
                    minRatio <= 0 || minRatio > 1) {
                stop("'minRatio' must be a single numeric value in [0,1]")
            }
            ratio[ratio < minRatio] <- minRatio
        }
    } else ratio <- NULL
    ## define local version of 'bwplot'
    localBwplot <- function(form, values, xlab = NULL, ylab = NULL, ..., 
            # these arguments are defined so that they aren't supplied twice:
            x, data, allow.multiple, outer, panel, groups) {
        bwplot(form, data=values, 
            panel=panelSpBwplot, xlab=xlab, ylab=ylab, ...)
    }
    ## call 'bwplot'
    localBwplot(as.formula(form), values, horizontal=horizontal, 
        coef=coef, zeros=zeros, ratio=ratio, do.out=FALSE, outliers=out, ...)
}


## panel function
panelSpBwplot <- function(x, y, coef=1.5, zeros = TRUE, 
        ratio, outliers, subscripts, ...) {
    out <- outliers[subscripts]
    if(zeros) {
        i <- packet.number()
        localPanelBwplot <- function(..., ratio, box.ratio, box.width) {
            panel.bwplot(..., box.ratio=ratio)
        }
        localPanelBwplot(x[!out], y[!out], coef=coef, ratio=ratio[,i], ...)
    } else panel.bwplot(x[!out], y[!out], coef=coef, ...)
    panel.points(x[out], y[out], ...)
}


## internal utility functions

# get data.frame and all required statistics
getBwplotStats <- function(x, weights = NULL, 
        cond = NULL, data, ..., name = "") {
    if(is.null(cond)) {
        x <- data[, x]
        w <- if(length(weights) == 0) NULL else data[, weights]
        prepBwplotStats(x, w, ..., name=name)
    } else {
        tmp <- tapply(1:nrow(data), data[, cond, drop=FALSE], 
            function(i) {
                x <- data[i, x]
                w <- if(length(weights) == 0) NULL else data[i, weights]
                g <- unique(data[i, cond, drop=FALSE])
                res <- prepBwplotStats(x, w, ..., name=name)
                res$values <- cbind(res$values, 
                    g[rep.int(1, nrow(res$values)), , drop=FALSE])
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
    values <- mapply(function(x, v) cbind(x$values, .var=v), 
        tmp, names(x), SIMPLIFY=FALSE, USE.NAMES=FALSE)
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

# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

spCdfplot <- function(x, ...) UseMethod("spCdfplot")

spCdfplot.default <- function(x, weights = NULL, cond = NULL, dataS, 
        dataP = NULL, approx = NULL, n = 10000, bounds = TRUE, ...) {
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
    nP <- if(inherits(dataP, "data.frame")) 2 else 1+length(dataP)
    # check 'approx'
    if(is.null(approx)) approx <- !is.null(dataP)
    else if(!is.logical(approx) || length(approx) == 0) approx <- FALSE
    if(is.null(dataP)) {
        if(length(approx) > 1) approx <- isTRUE(approx[1])
    } else {
        if(length(approx) == 1) {
            approx <- c(FALSE, rep.int(isTRUE(approx), nP-1))
        } else approx <- c(isTRUE(approx[1]), rep.int(isTRUE(approx[2]), nP-1))
    }
    # check 'n'
    if(any(approx) && (!is.numeric(n) || length(n) == 0)) {
        stop("'n' is not numeric or does not have positive length")
    }
    if(is.null(dataP)) {
        if(length(n) > 1) n <- if(approx) n[1] else NA
    } else {
        if(length(n) == 1) n <- ifelse(approx, n, NA)
        else n <- c(if(approx[1]) n[1] else NA, ifelse(approx[-1], n[2], NA))
    }
    # check 'bounds'
    bounds <- isTRUE(bounds)
    # define labels for grouping variable
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
    ## construct objects for 'xyplot'
    # from sample
    tmp <- getCdf(x, weights, cond, dataS, 
        approx=approx[1], n=n[1], name=lab[1])
    values <- tmp$values
    app <- t(tmp$approx)
    # from population(s)
    if(!is.null(dataP)) {
        if(inherits(dataP, "data.frame")) {
            tmp <- getCdf(x, NULL, cond, dataP, 
                approx=approx[2], n=n[2], name=lab[2])
            values <- rbind(values, tmp$values)
            app <- rbind(app, tmp$approx)
        } else {
            tmp <- mapply(function(dP, approx, n, l) {
                    getCdf(x, NULL, cond, dP, 
                        approx=approx, n=n, name=l)
                }, dataP, approx[-1], n[-1], lab[-1], 
                SIMPLIFY=FALSE, USE.NAMES=FALSE)
            values <- rbind(values, 
                do.call(rbind, lapply(tmp, function(x) x$values)))
            app <- rbind(app, do.call(rbind, lapply(tmp, function(x) x$approx)))
        }
    }
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
        # this produces a 'NOTE' during 'R CMD check':
#        xyplot(form, data=values, groups=if(nP == 1) NULL else .name, 
#            panel=panelSpCdfplot, prepanel=prepanelSpCdfplot, 
#            xlab=xlab, ylab=ylab, auto.key=auto.key, ...)
        command <- paste("xyplot(form, data=values,", 
            "groups=if(nP == 1) NULL else .name,", 
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
getCdf <- function(x, weights = NULL, 
        cond = NULL, data, ..., name = "") {
    if(is.null(cond)) {
        x <- data[, x]
        w <- if(length(weights) == 0) NULL else data[, weights]
        prepCdf(x, w, ..., name=name)
    } else {
        tmp <- tapply(1:nrow(data), data[, cond, drop=FALSE], 
            function(i) {
                x <- data[i, x]
                w <- if(length(weights) == 0) NULL else data[i, weights]
                g <- unique(data[i, cond, drop=FALSE])
                res <- prepCdf(x, w, ..., name=name)
                res$values <- cbind(res$values, 
                    g[rep.int(1, nrow(res$values)), , drop=FALSE])
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
    values <- mapply(function(x, v) cbind(x$values, .var=v), 
        tmp, names(x), SIMPLIFY=FALSE, USE.NAMES=FALSE)
    values <- do.call(rbind, values)
    approx <- sapply(tmp, function(x) x$approx)
    list(values=values, approx=approx)
}

prepCdf.default <- function(x, w, ..., name = "") {
    tmp <- spCdf(x, w, ...)
    values <- data.frame(.x=c(tmp$x[1], tmp$x), .y=c(0, tmp$y), .name=name)
    list(values=values, approx=tmp$approx)
}

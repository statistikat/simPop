# ----------------------------------------
# Authors: Andreas Alfons and Stefan Kraft
#          Vienna University of Technology
# ----------------------------------------

contingencyWt <- function(x, ...) UseMethod("contingencyWt")

contingencyWt.default <- function(x, y, weights = NULL, ...) {
    tab <- tableWt(data.frame(x, y), weights)
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0]
    chisq <- as.numeric(chisq.test(tab)$statistic)
    return(sqrt(chisq / (sum(tab) + chisq)))
}

#contingencyWt.default <- function(x, y, weights = NULL, ...) {
#    contingencyWt(data.frame(x, y), weights, ...)
#}

contingencyWt.matrix <- function(x, weights = NULL, ...) {
    contingencyWt(as.data.frame(x), weights=weights, ...)
}

contingencyWt.data.frame <- function(x, weights = NULL, ...) {
    # computes *pairwise* contingency coefficients
    p <- ncol(x)
    res <- matrix(NA, ncol=p-1, nrow=p-1)
    for(i in 1:(p-1)) {
        for(j in (i+1):p) {
            res[i, j-1] <- contingencyWt(x[, i], x[, j], weights)
        }
    }
    nam <- names(x)
    dimnames(res) <- list(nam[1:(p-1)], nam[2:p])
    res
}

#contingencyWt.data.frame <- function(x, weights = NULL, 
#        na.rm = FALSE, drop = FALSE, ...) {
#    # initializations
#    if(!is.data.frame(x)) x <- as.data.frame(x)
#    if(!is.null(weights)) {
#        if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
#        else if(length(weights) != nrow(x)) {
#            stop("length of 'weights' must equal the number of rows in 'x'")
#        } else if(!all(is.finite(weights))) stop("missing or infinite weights")
#    }
#    if(isTRUE(na.rm)) {
#        # apply is way too slow for large population data
##        ok <- !apply(is.na(x), 1, any)
#        # this is a fast but ugly replacement
#        # -----
#        command <- sapply(1:ncol(x), function(i) paste("is.na(x[, ", i, "])", sep=""))
#        command <- paste(command, collapse=" | ")
#        ok <- !eval(parse(text=command))
#        # -----
#        x <- x[ok, , drop=FALSE]
#        if(!is.null(weights)) weights <- weights[ok]
#    }
#    if(isTRUE(drop)) x <- dropLevels(x)
#    tab <- tableWt(x, weights)
#    chisq <- as.numeric(chisq.test(tab)$statistic)
#    return(sqrt(chisq / (sum(tab) + chisq)))
#}

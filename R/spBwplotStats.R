# ----------------------------------------
# Authors: Stefan Kraft and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

spBwplotStats <- function(x, weights = NULL, coef = 1.5, 
        zeros = TRUE, do.out = TRUE) {
    # initializations
    if(!is.numeric(x)) stop("'x' must be a numeric vector")
    if(!is.numeric(coef) || length(coef) != 1 || coef < 0) {
        stop("'coef' must be a single non-negative number")
    }
    # get quantiles
    if(isTRUE(zeros)) {
        zero <- ifelse(is.na(x), FALSE, x == 0)
        x <- x[!zero]
        if(is.null(weights)) nzero <- sum(zero)
        else {
            # if 'zeros' is not TRUE, these checks are done in 'quantileWt'
            # but here we need them since we use subscripting
            if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
            else if(length(weights) != length(zero)) {
                stop("'weights' must have the same length as 'x'")
            }
            nzero <- sum(weights[zero])
            weights <- weights[!zero]
        }
    } else nzero <- NULL
    ok <- !is.na(x)
    n <- if(is.null(weights)) sum(ok) else sum(weights[ok])
    if(n == 0) stats <- rep.int(NA, 5)
    else stats <- quantileWt(x, weights)
    iqr <- diff(stats[c(2, 4)])  # inter quartile range
    if(coef == 0) do.out <- FALSE
    else {
        if(is.na(iqr)) out <- is.infinite(x) 
        else {
            lower <- stats[2] - coef * iqr
            upper <- stats[4] + coef * iqr
            out <- ifelse(ok, x < lower | x > upper, FALSE)
        }
        if(any(out)) stats[c(1, 5)] <- range(x[!out], na.rm=TRUE)
    }
    res <- list(stats=stats, n=n, nzero=nzero, 
        out=if(isTRUE(do.out)) x[out] else numeric())
    class(res) <- "spBwplotStats"
    res
}

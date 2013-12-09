# ---------------------------------------
# Author: Stefan Kraft
#         Vienna University of Technology
# ---------------------------------------

quantileWt <- function(x, weights = NULL, 
        probs = seq(0, 1, 0.25), na.rm = TRUE) {
    # initializations
    if(!is.numeric(x)) stop("'x' must be a numeric vector")
    x <- unname(x)  # unlike 'quantile', this never returns a named vector
    if(is.null(weights)) {
        return(quantile(x, probs, na.rm=na.rm, names=FALSE, type=1))
    } else if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
    else if(length(weights) != length(x)) {
        stop("'weights' must have the same length as 'x'")
    } else if(!all(is.finite(weights))) stop("missing or infinite weights")
    if(!is.numeric(probs) || all(is.na(probs)) || 
            isTRUE(any(probs < 0 | probs > 1))) {
        stop("'probs' must be a numeric vector with values in [0,1]")
    }
    if(length(x) == 0) return(rep.int(NA, length(probs)))
    if(!isTRUE(na.rm) && any(is.na(x))) {
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
    return(q)
}

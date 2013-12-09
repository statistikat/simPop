# ----------------------------------------
# Authors: Stefan Kraft and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

spCdf <- function(x, weights = NULL, approx = FALSE, n = 10000) {
    ## initializations
    # remove non-finites values
    ok <- is.finite(x)
    x <- x[ok]
    m <- length(x)
    if(m == 0) stop("'x' must contain finite values")
    # order observations
    ord <- sort.list(x, na.last=NA, method = "quick")
    x <- x[ord]
    # check whether CDF should be approximated
    approx <- isTRUE(approx)
    if(approx) {
        if(!is.numeric(n) && length(n) != 1 && n <= 0) {
            stop("'n' must be a single positive integer")
        } else if(m < n) {
            approx <- FALSE
            warning("number of finite values in 'x' ", 
                "is smaller than 'n': no approximation")
        }
    }
    # define coordinates
    if(is.null(weights)) y <- (1:m)/m
    else {
        weights <- weights[ok][ord]
        cw <- cumsum(weights)
        y <- cw / cw[m]
    }
    res <- if(approx) approx(x, y, ties="ordered", n=n) else list(x=x, y=y)
    res$approx <- approx
    class(res) <- "spCdf"
    # return object
    res
}

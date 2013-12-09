# ----------------------------------------
# Authors: Andreas Alfons and Stefan Kraft
#          Vienna University of Technology
# ----------------------------------------

tableWt <- function(x, weights = NULL, useNA = c("no", "ifany", "always")) {
    # initializations
    if(!is.data.frame(x)) x <- as.data.frame(x)
    if(is.null(weights)) return(table(x, useNA=useNA))
    else if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
    else if(length(weights) != nrow(x)) {
        stop("length of 'weights' must equal the number of rows in 'x'")
    } else if(!all(is.finite(weights))) stop("missing or infinite weights")
    useNA <- match.arg(useNA)
    if(nrow(x) > 0 && ncol(x) > 0 && useNA != "no") {
        always <- useNA == "always"
        if(ncol(x) == 1) x[, 1] <- factorNA(x[, 1], always)
        else x <- as.data.frame(lapply(x, factorNA, always))
    }
    # compute and return weighted table
    tab <- round(tapply(weights, x, sum))
    tab[is.na(tab)] <- 0
    class(tab) <- "table"
    return(tab) 
}

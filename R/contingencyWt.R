#' Weighted contingency coefficients
#' 
#' Compute (weighted) pairwise contingency coefficients.
#' 
#' The function \code{\link{tableWt}} is used for the computation of the
#' corresponding pairwise contingency tables.
#' 
#' @name contingencyWt
#' @aliases contingencyWt contingencyWt.default contingencyWt.matrix
#' contingencyWt.data.frame
#' @docType method
#' @param x for the default method, a vector that can be interpreted as factor.
#' For the matrix and \code{data.frame} methods, the columns should be
#' interpretable as factors.
#' @param y a vector that can be interpreted as factor.
#' @param weights an optional numeric vector containing sample weights.
#' @param \dots for the generic function, arguments to be passed down to the
#' methods, otherwise ignored.
#' @return For the default method, the (weighted) contingency coefficient of
#' \code{x} and \code{y}.
#' 
#' For the matrix and \code{data.frame} method, a matrix of (weighted) pairwise
#' contingency coefficients for all combinations of columns.  Elements below
#' the diagonal are \code{NA}.
#' @author Andreas Alfons and Stefan Kraft
#' @export
#' @keywords methods
#' @seealso \code{\link{tableWt}}
#' @references Kendall, M.G. and Stuart, A. (1967) \emph{The Advanced Theory of
#' Statistics, Volume 2: Inference and Relationship}. Charles Griffin & Co Ltd,
#' London, 2nd edition.
#' @keywords category
#' @examples
#' 
#' data(eusilcS)
#' 
#' ## default method
#' contingencyWt(eusilcS$pl030, eusilcS$pb220a, weights = eusilcS$rb050)
#' 
#' ## data.frame method
#' basic <- c("age", "rb090", "hsize", "pl030", "pb220a")
#' contingencyWt(eusilcS[, basic], weights = eusilcS$rb050)
contingencyWt <- function(x, ...) UseMethod("contingencyWt")
NULL

#' @rdname contingencyWt
#' @export
contingencyWt.default <- function(x, y, weights = NULL, ...) {
    tab <- tableWt(data.frame(x, y), weights)
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0]
    chisq <- as.numeric(chisq.test(tab)$statistic)
    return(sqrt(chisq / (sum(tab) + chisq)))
}
NULL

#contingencyWt.default <- function(x, y, weights = NULL, ...) {
#    contingencyWt(data.frame(x, y), weights, ...)
#}

#' @rdname contingencyWt
#' @export
contingencyWt.matrix <- function(x, weights = NULL, ...) {
    contingencyWt(as.data.frame(x), weights=weights, ...)
}
NULL

#' @rdname contingencyWt
#' @export
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
NULL

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

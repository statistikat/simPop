# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

spTable <- function(dataS, dataP, select = NULL, weights = NULL) {
    UseMethod("spTable", select)
}

spTable.formula <- function(dataS, dataP, select, weights = NULL) {
    # prepare formula
    tmp <- as.character(select)
    tmp <- tmp[length(tmp)]  # omit LHS
    tmp <- sub("\\|.*", "", tmp)  # RHS (omit conditioning variables)
    f <- formula(c("~", tmp))  # formula
    if(any(attr(terms(f, data=dataS), "order") > 1)) {
        stop("interactions are not allowed in the formula")
    }
    weights <- eval(substitute(weights), dataS, environment(select))  # evaluate
    # data sets
    dataS <- get_all_vars(f, dataS)  # does not work correctly 
    dataP <- get_all_vars(f, dataP)  # if, e.g., I(x+y) is used
    # call default method
    spTable(dataS, dataP, weights=weights)
}

spTable.default <- function(dataS, dataP, select = NULL, weights = NULL) {
    # initializations
    if(!is.data.frame(dataS)) stop("'dataS' must be a data.frame")
    if(!is.data.frame(dataP)) stop("'dataP' must be a data.frame")
    # prepare weights
    if(is.null(weights)) {
        n <- nrow(dataS)
        weights <- rep.int(nrow(dataP)/n, n)
    } else if(is.character(weights) && length(weights) == 1) {
        weights <- dataS[, weights]
        if(is.null(select)) dataS <- dataS[, setdiff(names(dataS), weights)]
    } else if(!is.numeric(weights)) {
        stop("'weights' must be a character string ", 
            "specifying a column of 'dataS' or a numeric vector")
    }
    # prepare data seta
    if(!is.null(select)) {
        if(!is.character(select)) {
            stop("'select' must be a character vector ", 
                "specifying columns of 'dataS' and 'dataP'")
        }
        dataS <- dataS[, select, drop=FALSE]
        dataP <- dataP[, select, drop=FALSE]
    }
    # compute weighted table from sample (expected)
    tabS <- tableWt(dataS, weights)
    # compute table for population (realized)
    tabP <- table(dataP)
    # create object and return result
    res <- list(expected=tabS, realized=tabP)
    class(res) <- "spTable"
    res
}


# methods for class "spTable"

as.array.spTable <- function(x, ...) {
    values <- c(as.integer(x$expected), as.integer(x$realized))
    d <- c(dim(x$expected), 2)
    dnNew <- list(which=c("expected", "realized"))
#    dnNew <- list(Data=c("Sample", "Population"))
    dn <- c(dimnames(x$expected), dnNew)
    array(values, dim=d, dimnames=dn)
}

as.table.spTable <- function(x, ...) {
    tab <- as.array(x)
    class(tab) <- "table"
    tab
}

plot.spTable <- function(x, ...) spMosaic(x, ...)

print.spTable <- function(x, ...) {
    # expected (from sample)
    cat("Expected:\n")
    print(x$expected, ...)
    # realized (population)
    cat("\nRealized:\n")
    print(x$realized, ...)
}

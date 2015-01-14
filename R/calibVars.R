calibVars <- function(x) UseMethod("calibVars")
calibVars.default <- function(x) {
  if(length(x) == 0) matrix(integer(), 0, 0)
  x <- as.factor(x)
  levs <- levels(x)
  res <- .Call("simPop_binary_representation", levels=1:length(levels(x)), values=as.integer(x), package="simPop")
  colnames(res) <- levs
  rownames(res) <- names(x)  # set rownames from original vector
  res    
}
calibVars.matrix <- function(x) calibVars(as.data.frame(x))
calibVars.data.frame <- function(x) {
  res <- lapply(x, calibVars)  # list of matrices for each variable
  res <- mapply(function(x, nam) {
    colnames(x) <- paste(nam, colnames(x), sep=".")
    x
  }, res, names(x), SIMPLIFY=FALSE)
  res <- do.call("cbind", res)
  rownames(res) <- row.names(x)
  res
}

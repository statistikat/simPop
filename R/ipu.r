ipu <- function(inp, con, hid=NULL, eps=1e-07, verbose=FALSE) {
  if ( class(inp)[1] %in% c("data.frame", "data.table") ) {
    cnames <- colnames(inp)
    inp <- as.matrix(inp)
    colnames(inp) <- cnames
  }

  if ( is.null(hid) ) {
    hhid <- 1:nrow(inp)
  } else {
    ii <- match(hid, colnames(inp))
    hhid <- inp[,ii]
    inp <- inp[,-c(ii)]
  }
  inp <- inp[,match(names(con), colnames(inp))]
  con <- con[match(colnames(inp), names(con))]
  if ( any(colnames(inp) != names(con)) ) {
    stop("constraints (con) do not match input (inp) -> check your input!\n")
  }
  w <- as.numeric(rep(1,nrow(inp)))
  con <- as.numeric(unlist(con))
  w <- .Call("simPop_ipu_work", inp, con, w, eps=eps, verbose=ifelse(verbose, 1L, 0L), package="simPop")
  out <- cbind(hhid, inp)
  out <- as.data.frame(out)
  out$weights <- w
  invisible(out)
}

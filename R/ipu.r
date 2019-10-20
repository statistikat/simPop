#' iterative proportional updating
#'
#' adjust sampling weights to given totals based on household-level and/or
#' individual level constraints
#'
#' @name ipu
#' @param inp a \code{data.frame} or \code{data.table} containing household ids
#' (optionally), counts for household and/or personal level attributes that
#' should be fitted.
#' @param con named list with each list element holding a constraint total with
#' list-names relating to column-names in \code{inp}.
#' @param hid character vector specifying the variable containing household-ids
#' within \code{inp} or NULL if such a variable does not exist.
#' @param eps number specifiying convergence limit
#' @param verbose if TRUE, ipu will print some progress information.
#' @export
#' @author Bernhard Meindl
#' @keywords method
#' @examples
#' # basic example
#' inp <- as.data.frame(matrix(0, nrow=8, ncol=6))
#' colnames(inp) <- c("hhid","hh1","hh2","p1","p2","p3")
#' inp$hhid <- 1:8
#' inp$hh1[1:3] <- 1
#' inp$hh2[4:8] <- 1
#' inp$p1 <- c(1,1,2,1,0,1,2,1)
#' inp$p2 <- c(1,0,1,0,2,1,1,1)
#' inp$p3 <- c(1,1,0,2,1,0,2,0)
#' con <- list(hh1=35, hh2=65, p1=91, p2=65, p3=104)
#' res <- ipu(inp=inp, hid="hhid", con=con, verbose=FALSE)
#'
#' # more sophisticated
#' # load sample and population data
#' data(eusilcS)
#' data(eusilcP)
#'
#' # variable generation and preparation
#' eusilcS$hsize <- factor(eusilcS$hsize)
#'
#' # make sure, factor levels in sample and population match
#' eusilcP$region <- factor(eusilcP$region, levels = levels(eusilcS$db040))
#' eusilcP$gender <- factor(eusilcP$gender, levels = levels(eusilcS$rb090))
#' eusilcP$hsize  <- factor(eusilcP$hsize , levels = levels(eusilcS$hsize))
#'
#' # generate input matrix
#' # we want to adjust to variable "db040" (region) as household variables and
#' # variable "rb090" (gender) as individual information
#' library(data.table)
#' samp <- data.table(eusilcS)
#' pop <-  data.table(eusilcP)
#' setkeyv(samp, "db030")
#' hh <- samp[!duplicated(samp$db030),]
#' hhpop <- pop[!duplicated(pop$hid),]
#'
#' # reg contains for each region the number of households
#' reg <- data.table(model.matrix(~db040 +0, data=hh))
#' # hsize contains for each household size the number of households
#' hsize <- data.table(model.matrix(~factor(hsize) +0, data=hh))
#'
#' # aggregate persons-level characteristics per household
#' # gender contains for each household the number of males and females
#' gender <- data.table(model.matrix(~db030+rb090 +0, data=samp))
#' setkeyv(gender, "db030")
#' gender <- gender[, lapply(.SD, sum), by = key(gender)]
#'
#' # bind together and use it as input
#' inp <- cbind(reg, hsize, gender)
#'
#' # the totals we want to calibrate to
#' con <- c(
#'   as.list(xtabs(rep(1, nrow(hhpop)) ~ hhpop$region)),
#'   as.list(xtabs(rep(1, nrow(hhpop)) ~ hhpop$hsize)),
#'   as.list(xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$gender))
#' )
#' # we need to have the same names as in 'inp'
#' names(con) <- setdiff(names(inp), "db030")
#'
#' # run ipu und check results
#' res <- ipu(inp=inp, hid="db030", con=con, verbose=TRUE)
#'
#' is <- sapply(2:(ncol(res)-1), function(x) {
#'   sum(res[,x]*res$weights)
#' })
#' data.frame(required=unlist(con), is=is)
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
  w <- ipu_work(
    inp = inp,
    con = con,
    w = w,
    eps = eps,
    verbose = ifelse(verbose, 1L, 0L)
  )
  out <- cbind(hhid, inp)
  out <- as.data.frame(out)
  out$weights <- w
  invisible(out)
}

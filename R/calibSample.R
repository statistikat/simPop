calibSample <- function(X, d, totals, q=NULL,
  method=c("raking", "linear", "logit"),
  bounds=c(0, 10), maxit=500, tol=1e-06,
  eps=.Machine$double.eps) {

  ## initializations and error handling
  X <- as.matrix(X)
  d <- as.numeric(d)
  totals <- as.numeric(totals)
  haveNA <- c(any(is.na(X)), any(is.na(d)),
    any(is.na(totals)), !is.null(q) && any(is.na(q)))
  if ( any(haveNA) ) {
    argsNA <- c("'X'", "'d'", "'totals'", "'q'")[haveNA]
    stop("missing values in the following arguments", paste(argsNA, collapse=", "),"\n")
  }
  n <- nrow(X)  # number of rows
  if ( length(d) != n ) {
    stop("length of 'd' not equal to number of rows in 'X'!\n")
  }
  p <- ncol(X)  # number of columns
  if ( length(totals) != p ) {
    stop("length of 'totals' not equal to number of columns in 'X'!\n")
  }
  if ( is.null(q) ) {
    q <- rep.int(1, n)
  } else {
    q <- as.numeric(q)
    if ( length(q) != n ) {
      stop("length of 'q' not equal to number of rows in 'X'!\n")
    }
    if ( any(is.infinite(q)) ) {
      stop("infinite values in 'q'")
    }
  }
  method <- match.arg(method)

  ## computation of g-weights
  if ( method == "linear" ) {
    ## linear method (no iteration!)
    lambda <- ginv(t(X * d * q) %*% X, tol=eps) %*% (totals - as.vector(t(d) %*% X))
    g <- 1 + q * as.vector(X %*% lambda)  # g-weights
  } else {
    ## multiplicative method (raking) or logit method
    lambda <- matrix(0, nrow=p)  # initial values
    # function to determine whether teh desired accuracy has
    # not yet been reached (to be used in the 'while' loop)
    tolNotReached <- function(X, w, totals, tol) {
      max(abs(crossprod(X, w) - totals)/totals) >= tol
    }
    if ( method == "raking" ) {
      ## multiplicative method (raking)
      # some initial values
      g <- rep.int(1, n)  # g-weights
      w <- d  # sample weights
      ## iterations
      i <- 1
      while ( !any(is.na(g)) && tolNotReached(X, w, totals, tol) && i <= maxit ) {
        # here 'phi' describes more than the phi function in Deville,
        # Saerndal and Sautory (1993); it is the whole last term of
        # equation (11.1)
        phi <- t(X) %*% w - totals
        T <- t(X * w)
        dphi <- T %*% X  # derivative of phi function (to be inverted)
        lambda <- lambda - ginv(dphi, tol=eps) %*% phi  # update 'lambda'
        g <- exp(as.vector(X %*% lambda) * q)  # update g-weights
        w <- g * d  # update sample weights
        i <- i + 1  # increase iterator
      }
      ## check wether procedure converged
      if( any(is.na(g)) || i > maxit ) {
        warning("no convergence!\n")
        g <- NULL
      }
    } else {
      ## logit (L, U) method
      ## error handling for bounds
      if ( length(bounds) < 2 ) {
        stop("'bounds' must be a vector of length 2")
      } else {
        bounds <- bounds[1:2]
      }
      if ( bounds[1] >= 1 ) {
        stop("the lower bound must be smaller than 1")
      }
      if ( bounds[2] <= 1 ) {
        stop("the lower bound must be larger than 1")
      }
      ## some preparations
      A <- diff(bounds)/((1 - bounds[1]) * (bounds[2] - 1))
      # function to bound g-weights
      getG <- function(u, bounds) {
        (bounds[1] * (bounds[2]-1) + bounds[2] * (1-bounds[1]) * u) /
        (bounds[2]-1 + (1-bounds[1]) * u)
      }
      ## some initial values
      g <- getG(rep.int(1, n), bounds)  # g-weights
      # in the procedure, g-weights outside the bounds are moved to the
      # bounds and only the g-weights within the bounds are adjusted.
      # these duplicates are needed since in general they are changed in
      # each iteration while the original values are also needed
      X1 <- X
      d1 <- d
      totals1 <- totals
      q1 <- q
      g1 <- g
      indices <- 1:n
      # function to determine which g-weights are outside the bounds
      anyOutOfBounds <- function(g, bounds) {
        any(g < bounds[1]) || any(g > bounds[2])
      }
      ## iterations
      i <- 1
      while ( !any(is.na(g)) && (tolNotReached(X, g*d, totals, tol ) || anyOutOfBounds(g, bounds)) && i <= maxit) {
        # if some of the g-weights are outside the bounds, these values
        # are moved to the bounds and only the g-weights within the
        # bounds are adjusted
        if ( anyOutOfBounds(g, bounds) ) {
          g[g < bounds[1]] <- bounds[1]
          g[g > bounds[2]] <- bounds[2]
          # values within the bounds
          tmp <- which(g > bounds[1] & g < bounds[2])
          if ( length(tmp) > 0 ) {
            indices <- tmp
            X1 <- X[indices,]
            d1 <- d[indices]
            if ( length(indices) < n ) {
              totals1 <- totals - as.vector(t(g[-indices] * d[-indices]) %*% X[-indices, , drop=FALSE])
            }
            q1 <- q[indices]
            g1 <- g[indices]
          }
        }
        w1 <- g1 * d1  # current sample weights
        # here 'phi' describes more than the phi function in Deville,
        # Saerndal and Sautory (1993); it is the whole last term of
        # equation (11.1)
        phi <- t(X1) %*% w1 - totals1
        T <- t(X1 * w1)
        dphi <- T %*% X1  # derivative of phi function (to be inverted)
        lambda <- lambda - ginv(dphi, tol=eps) %*% phi  # update 'lambda'
        # update g-weights
        u <- exp(A * as.vector(X1 %*% lambda) * q1)
        g1 <- getG(u, bounds)
        g[indices] <- g1
        i <- i+1  # increase iterator
      }
      ## check wether procedure converged
      if ( any(is.na(g)) || i > maxit ) {
        warning("no convergence!\n")
        g <- NULL
      }
    }
  }
  ## return g-weights
  invisible(g)
}

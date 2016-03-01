#' Correct 5-year age heaping
#' 
#' Correct for age heaping using truncated (log-)normal distributions
#' 
#' @details Age heaping can cause substantial bias
#' in important measures and thus age heaping should be corrected. 
#' 
#' For method \dQuote{lnorm}, a truncated log-normal is fit to the whole age distribution. 
#' Then for each age heap (at 0, 5, 10, 15, ...) random numbers of a truncated 
#' log-normal (with lower and upper bound) is drawn in the 
#' interval +- 2 around the heap (rounding of degree 2) 
#' using the inverse transformation method. A ratio of randomly 
#' chosen observations on an age heap are replaced 
#' by these random draws. For the ratio the age distribution is chosen, whereas 
#' on an age heap (e.g. 5)
#' the arithmetic means of the two neighboring ages are calculated 
#' (average counts on age 4 and age 6 for age heap equals 5, for example).
#' The ratio on, e.g. age equals 5 is then given by the count on age 5 divided by this mean
#' This is done for any age heap at (0, 5, 10, 15, ...). 
#' 
#' Method \dQuote{norm} replace the draws from truncated log-normals to draws from 
#' truncated normals. It depends on the age distrubution (if right-skewed or not) if method
#' \dQuote{lnorm} or \dQuote{norm} should be used. Many distributions with heaping problems
#' are right-skewed.
#' 
#' Method \dQuote{unif} draws the mentioned ratio of observations on truncated uniform distributions
#' around the age heaps.
#' 
#' Repeated calls of this function mimics multiple imputation, i.e. repeating this 
#' procedure m times provides m imputed datasets that properly reflect the 
#' uncertainty from imputation.
#' @param x numeric vector
#' @param method default is log-normal (\dQuote{lnorm}). Method \dQuote{norm} and \dQuote{unif} are
#' also available.
#' @author Matthias Templ
#' @export
#' @examples 
#' ## test with toy data (to see which methods works better)
#' ## meanlog and sdlog estimated from real data.
#' ## age distribution without heaping problems:
#' y <- rlnorm(10000, meanlog = 2.466869, sdlog = 1.652772)
#' y <- round(y[y < 90])
#' plot(table(y))
#' ## artificially introduce age heaping:
#' w <- y %in% seq(0,100, by=5)
#' z <- c(y, y[w])
#' z <- sample(z)
#' ## same size as y:
#' z <- z[1:length(y)]
#' plot(table(z))
#' ## correct age heaping:
#' a1 <- correctHeap(z)
#' a2 <- correctHeap(z, method = "norm")
#' a3 <- correctHeap(z, method = "unif")
#' plot(table(a1))
#' plot(table(a2))
#' plot(table(a3))
#' sum(abs(table(a1) - table(y)))
#' sum(abs(table(a2) - table(y)))
#' sum(abs(table(a3) - table(y)))
correctHeap <- function(x, method = "lnorm"){
  s <- seq(0, 95, 5) 
  w <- x %in% s
  tab <- table(x)
  keep <- sapply(s + 1, function(x) mean(c(tab[x-1], tab[x+1])))
  ratio <- tab[s + 1] / keep
  if( method == "lnorm"){
    age0 <- as.numeric(x)
    age0[age0 == 0] <- 0.01
    logn <- fitdistrplus::fitdist(age0, "lnorm")
    j <- 0
    for(i in s[-length(s)]){
      j <- j + 1
      index <- which(x == i)
      if(!is.na(ratio[j])){
        r <- sample(index, size = max(c(1, length(index) - length(index) / ratio[j])))
        x[r] <- round(EnvStats::rlnormTrunc(length(r), 
                                  meanlog = logn$estimate[1], 
                                  sdlog = as.numeric(logn$estimate[2]), 
                                  min = max(i-2, 0), max = i+2))
      }
    }
  }
  if( method == "norm"){
    j <- 0
    for(i in s[-length(s)]){
      j <- j + 1
      index <- which(x == i)
      if(!is.na(ratio[j])){
        r <- sample(index, size = max(c(1, length(index) - length(index) / ratio[j])))
        x[r] <- round(EnvStats::rnormTrunc(length(r), 
                                 mean = i, 
                                 sd = 1, 
                                 min = max(i-2, 0), max = i+2))
      }
    }    
  }
  if( method == "unif"){
    j <- 0
    for(i in s[-length(s)]){
      j <- j + 1
      index <- which(x == i)
      if(!is.na(ratio[j])){
        r <- sample(index, size = max(c(1, length(index) - length(index) / ratio[j])))
        x[r] <- sample(max(c(i-2,0)):(i+2), length(r), replace = TRUE)
      }
    }      
  }
  x
}




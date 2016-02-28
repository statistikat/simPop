#' Correct 5-year age heaping
#' 
#' Correct fo age heaping using truncated (log-)normal distributions
#' 
#' @details for method lnorm, a truncated log-normal is fit to the whole age distribution. Then
#' for each age heaps (0,5,10,15, ...) random numbers of truncated log-normal is drawn in the 
#' interval +- 2 of the age heap for a certain ratio of values at an age heap. 
#' This ratio is taken from the mean of +- 1 of age regarding to the value at the age heap.
#' @param x numeric vector
#' @param method default is log-normal (\dQuote{lnorm}). Method \dQuote{norm} and \dQuote{unif} are
#' also available.
#' @author Matthias Templ
#' @export
#' @examples 
#' ## test with toy data (to see which methods works better)
#' ## meanlog and sdlog estimated from real data.
#' y <- rlnorm(10000, meanlog = 2.466869, sdlog = 1.652772)
#' y <- round(y[y < 90])
#' plot(table(y))
#' w <- y %in% seq(0,100, by=5)
#' z <- c(y, y[w])
#' z <- sample(z)
#' ## same size as y
#' z <- z[1:length(y)]
#' plot(table(z))
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




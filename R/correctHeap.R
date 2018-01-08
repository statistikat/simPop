#' Correct  age heaping
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
#' @param heaps
#' \itemize{
#' \item \code{5year}: heaps are assumed to be every 5 years (0,5,10,...)
#' \item \code{10year}: heaps are assumed to be every 10 years (0,10,20,...)
#' }
#' @param method a character specifying the algorithm used to correct the age
#' heaps. Allowed values are
#' \itemize{
#' \item \code{lnorm}: drawing from a truncated log-normal distribution. The
#' required parameters are estimated using original input data.
#' \item \code{norm}: drawing from a truncated normal distribution. The
#' required parameters are estimated using original input data.
#' \item \code{unif}: random sampling from a (truncated) uniform distribution
#' }
#' @param start a numeric value for the starting of the 5 or 10 year sequences
#' (e.g. 0, 5 or 10)
#' @param fixed numeric index vector with observation that should not be changed
#' @author Matthias Templ, Bernhard Meindl, Alexander Kowarik
#' @references 
#' M. Templ, B. Meindl, A. Kowarik, A. Alfons, O. Dupriez (2017) Simulation of Synthetic Populations for Survey Data Considering Auxiliary
#' Information. \emph{Journal of Statistical Survey}, \strong{79} (10), 1--38. doi: 10.18637/jss.v079.i10
#' @examples 
#' ## create some artificial data
#' age <- rlnorm(10000, meanlog=2.466869, sdlog=1.652772)
#' age <- round(age[age < 93])
#' barplot(table(age))
#'
#' ## artificially introduce age heaping and correct it:
#' # heaps every 5 years
#' year5 <- seq(0, max(age), 5)
#' age5 <- sample(c(age, age[age %in% year5]))
#' cc5 <- rep("darkgrey", length(unique(age)))
#' cc5[year5+1] <- "yellow"
#' barplot(table(age5), col=cc5)
#' barplot(table(correctHeaps(age5, heaps="5year", method="lnorm")), col=cc5)
#'
#' # heaps every 10 years
#' year10 <- seq(0, max(age), 10)
#' age10 <- sample(c(age, age[age %in% year10]))
#' cc10 <- rep("darkgrey", length(unique(age)))
#' cc10[year10+1] <- "yellow"
#' barplot(table(age10), col=cc10)
#' barplot(table(correctHeaps(age10, heaps="10year", method="lnorm")), col=cc10)
#' 
#' # the first 5 observations should be unchanged
#' barplot(table(correctHeaps(age10, heaps="10year", method="lnorm", fixed=1:5)), col=cc10)
#'

#' @return a numeric vector without age heaps
#' @export
correctHeaps <- function(x, heaps="10year", method="lnorm",start=0, fixed=NULL) {
  if ( !method %in% c("lnorm","norm","unif")) {
    stop("unsupported value in argument 'method'!\n")
  }
  if ( !heaps %in% c("5year","10year") ) {
    stop("unsupported value in argument 'heaps'!\n")
  }
  if(start>max(x)){
    stop("Starting Ageyear is greater than the maximum age in the data.")
  }
  if ( heaps=="10year" ) {
    s <- seq(start, max(x), by=10)
  }
  if ( heaps=="5year" ) {
    s <- seq(start, max(x), by=5)
  }
  
  tab <- table(x)
  keep <- sapply(s+1, function(x) mean(c(tab[x-1], tab[x+1])))
  ratio <- tab[s+1] / keep
  
  if ( method=="lnorm" ) {
    age0 <- as.numeric(x)
    age0[age0 == 0] <- 0.01
    logn <- fitdist(age0, "lnorm")
  }
  
  for ( j in 1:length(s) ) {
    i <- s[j]
    index <- which(x == i)
    if ( is.na(ratio[j]) ) {
      ssize <- 0
    } else {
      ssize <- ceiling(length(index) - length(index) / ratio[j])
    }
    if ( ssize>0 ) {
      # we need to take care of odd-years between leaps
      # thus we need to sample twice
      if ( heaps=="10year" ) {
        size1 <- ceiling(ssize/2)
        if(is.null(fixed)){
          r1 <- sample(index, size=size1)
        }else{
          indexTmp <- index[!index%in%fixed]
          if(length(indexTmp)==0){
            warning("There is no suitable observation to be changed left.")
            r1 <- NULL
          }else{
            r1 <- sample(indexTmp, size=min(c(length(indexTmp),size1)))  
          }
        }
        
        
        llow <- max(c(i-4,0))
        lup <- min((i+4),max(x))
        if(length(r1)>0){
          if ( method=="lnorm") {
            x[r1] <- round(rlnormTrunc(length(r1),
                                       meanlog=logn$estimate[1],
                                       sdlog=as.numeric(logn$estimate[2]),
                                       min=llow, max=lup))
          }
          if ( method=="norm") {
            x[r1] <- round(rnormTrunc(length(r1),
                                      mean=i, sd=1, min=llow, max=lup))
          }
          if ( method=="unif") {
            x[r1] <- sample(llow:lup, length(r1), replace=TRUE)
          }
          # sample 2:
          size2 <- ssize-size1
          if(is.null(fixed)){
            r2 <- sample(setdiff(index,r1), size=size2)  
          }else{
            indexTmp2 <- setdiff(indexTmp,r1)
            if(length(indexTmp2)==0){
              warning("There is no suitable observation to be changed left.")
              r2 <- NULL
            }else{
              r2 <- sample(indexTmp2, size=min(c(length(indexTmp2),size2)))  
            }
            
          }
          
          llow <- max(c(i-5,0))
          lup <- min((i+5),max(x))
          if(length(r2)>0){
            if ( method=="lnorm") {
              x[r2] <- round(rlnormTrunc(length(r2),
                                         meanlog=logn$estimate[1],
                                         sdlog=as.numeric(logn$estimate[2]),
                                         min=llow, max=lup))
            }
            if ( method=="norm") {
              x[r2] <- round(rnormTrunc(length(r2),
                                        mean=i, sd=1, min=llow, max=lup))
            }
            if ( method=="unif") {
              x[r2] <- sample(llow:lup, length(r2), replace=TRUE)
            }
          }
          
        }
      }
      if ( heaps=="5year" ) {
        if(is.null(fixed)){
          r <- sample(index, size=ssize)  
        }else{
          indexTmp <- index[!index%in%fixed]
          if(length(indexTmp)==0){
            warning("There is no suitable observation to be changed left.")
            r <- NULL
          }else{
            r <- sample(indexTmp, size=min(c(length(indexTmp),ssize)))  
          } 
        }
        
        llow <- max(c(i-2,0))
        lup <- min((i+2),max(x))
        if(length(r)>0){
          if ( method=="lnorm" ) {
            x[r] <- round(rlnormTrunc(length(r),
                    meanlog=logn$estimate[1],
                    sdlog=as.numeric(logn$estimate[2]),
                    min=llow, max=lup))
          }
          if ( method=="norm" ) {
            x[r] <- round(rnormTrunc(length(r),
                    mean=i, sd=1, min=llow, max=lup))
          }
          if ( method=="unif" ) {
            x[r] <- sample(llow:lup, length(r), replace=TRUE)
          }
        }
      }
    }
  }
  return(x)
}

#' correctSingleHeap
#'
#' Correct a specific age heap in a vector containing age in years
#'
#' @param x numeric vector representing age in years (integers)
#' @param heap numeric or integer vector of length 1 specifying the year
#' for which a heap should be corrected
#' @param before numeric or integer vector of length 1 specifying the number
#' of years before the heap that may be used to correct the heap. This input will
#' be rounded!
#' @param after numeric or integer vector of length 1 specifying the number
#' of years after the heap that may be used to correct the heap. This input will
#' be rounded!
#' \itemize{
#' \item \code{5year}: heaps are assumed to be every 5 years (0,5,10,...)
#' \item \code{10year}: heaps are assumed to be every 10 years (0,10,20,...)
#' }
#' @param method a character specifying the algorithm used to correct the age
#' heaps. Allowed values are
#' \itemize{
#' \item \code{lnorm}: drawing from a truncated log-normal distribution. The
#' required parameters are estimated using original input data.
#' \item \code{norm}: drawing from a truncated normal distribution. The
#' required parameters are estimated using original input data.
#' \item \code{unif}: random sampling from a (truncated) uniform distribution
#' }
#' @param fixed numeric index vector with observation that should not be changed
#' 
#' @author Matthias Templ, Bernhard Meindl, Alexander Kowarik
#' @return a numeric vector without age heaps
#' @export
#' @examples
#' ## create some artificial data
#' age <- rlnorm(10000, meanlog=2.466869, sdlog=1.652772)
#' age <- round(age[age < 93])
#' barplot(table(age))
#'
#' ## artificially introduce an age heap for a specific year
#' ## and correct it
#' age23 <- c(age, rep(23, length=sum(age==23)))
#' cc23 <- rep("darkgrey", length(unique(age)))
#' cc23[24] <- "yellow"
#' barplot(table(age23), col=cc23)
#' barplot(table(correctSingleHeap(age23, heap=23, before=2, after=3, method="lnorm")), col=cc23)
#' barplot(table(correctSingleHeap(age23, heap=23, before=5, after=5, method="lnorm")), col=cc23)
#' 
#' # the first 5 observations should be unchanged
#' barplot(table(correctSingleHeap(age23, heap=23, before=5, after=5, method="lnorm",
#'   fixed=1:5)), col=cc23)
correctSingleHeap <- function(x, heap, before=2, after=2, method="lnorm", fixed=NULL) {
  i <- NULL
  if ( !method %in% c("lnorm","norm","unif")) {
    i
    stop("unsupported value in argument 'method'!\n")
  }
  if ( !heap %in% unique(x) ) {
    stop("specified heap is not available in 'x'.\n")
  }
  if ( length(heap)!=1 ) {
    stop("you can specify only one heap!\n")
  }
  if ( !is.numeric(before) | !is.numeric(after) ) {
    stop("arguments 'before' and 'after' must be numeric!")
  }
  if ( length(before)!=1 | length(after)!=1 ) {
    stop("arguments 'before' and 'after' must be vectors with one element each!\n")
  }
  if ( before<0 | after < 0 ) {
    stop("arguments 'before' and 'after' must be positive numbers!\n")
  }
  before <- round(before)
  after <- round(after)
  
  llow <- heap-before
  lup <- heap+after
  if ( llow < 0 | lup > max(x) ) {
    stop("paramters 'before' or 'after' are too large!\n")
  }
  
  tab <- table(x)
  keep <- sapply(heap+1, function(x) mean(c(tab[x-1], tab[x+1])))
  ratio <- tab[heap+1] / keep
  
  if ( method=="lnorm" ) {
    age0 <- as.numeric(x)
    age0[age0 == 0] <- 0.01
    logn <- fitdist(age0, "lnorm")
  }
  
  index <- which(x == heap)
  if ( is.na(ratio) ) {
    ssize <- 0
  } else {
    ssize <- ceiling(length(index) - length(index) / ratio)
  }
  if ( ssize>0 ) {
    if(is.null(fixed)){
      r <- sample(index, size=ssize)
    }else{
      indexTmp <- index[!index%in%fixed]
      if(length(indexTmp)==0){
        warning("There is no suitable observation to be changed left.")
        r <- NULL
      }else{
        r <- sample(indexTmp, size=min(c(length(indexTmp),ssize)))  
      }
      
    }
    if(length(r)>0){
      if ( method=="lnorm" ) {
        x[r] <- round(rlnormTrunc(length(r),
                                  meanlog=logn$estimate[1],
                                  sdlog=as.numeric(logn$estimate[2]),
                                  min=llow, max=lup))
      }
      if ( method=="norm" ) {
        x[r] <- round(rnormTrunc(length(r),
                                 mean=heap, sd=1, min=llow, max=lup))
      }
      if ( method=="unif" ) {
        x[r] <- sample(llow:lup, length(r), replace=TRUE)
      }
    }
  }
  return(x)
}

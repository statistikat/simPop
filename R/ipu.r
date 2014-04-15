ipu <- function(inp, con, hid=NULL, eps=1e-07, verbose=FALSE) {
  # check for existing hid and possibly generate it
  # attach hhid later  
  if ( is.null(hid) ) {
    hhid <- 1:nrow(inp)
  } else {
    ii <- match(hid, colnames(inp))
    hhid <- inp[,ii,with=F]
    inp <- inp[,-c(ii), with=FALSE]      
  }
  
  # remove cols from inp not existing in con
  inp <- inp[,match(names(con), colnames(inp)),with=F]
  
  # reorder constraints
  con <- con[match(colnames(inp), names(con))]  
  if ( any(colnames(inp) != names(con)) ) {
    stop("constraints (con) does not match input (inp) -> check your input!\n")
  }
  
  # initialize weights
  w <- rep(1,nrow(inp))  
  gamma.vals <- sapply(1:length(con), function(x) { (abs(sum(w*inp[,x]) - con[[x]])) / con[[x]] })
  gamma <- mean(gamma.vals)
  
  run_ind <- TRUE
  counter <- 0
  while ( run_ind ) {
    counter <- counter + 1
    for ( j in 1:length(con) ) {
      v <- con[[j]] / sum(inp[,j, with=FALSE]*w, na.rm=TRUE)  
      index <- which(inp[,j, with=F]!=0)
      w[index] <- v*w[index]      
    }    
    gamma.vals.new <- sapply(1:length(con), function(x) { (abs(sum(w*inp[,x,with=F]) - con[[x]])) / con[[x]] })
    gamma.new <- mean(gamma.vals.new)
    delta <- abs(gamma.new - gamma)
    if ( verbose ) {
      cat("improvement in run",counter,":",delta," | gamma=",gamma.new,"\n")
    }
    
    if ( gamma.new < eps || delta < eps ) {
      run_ind <- FALSE
    } else {
      gamma.vals <- gamma.vals.new
      gamma <- gamma.new
    }
  }
  if ( verbose ) {
    cat("ipu finished after",counter,"interations!\n")
  }  
  out <- cbind(hhid, inp)  
  out$weights <- w
  invisible(out)
}

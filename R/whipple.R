whipple <- function(x, method="standard"){
  x <- x[x >= 23 & x <= 62]
  if(method == "standard"){
    xm <- x %% 5
    whipple <- (length(xm[xm==0])/length(x))*500
  }
  if(method == "modified"){
    tab <- table(x)
    sp <- function(p) seq(p,p+30,10)
    W <- numeric(9)
    W[1] <- 5*sum(sp(31)) / sum(5*sp(29)) 
    W[2] <- 5*sum(sp(32)) / sum(5*sp(30)) 
    W[3] <- 5*sum(sp(23)) / sum(5*sp(21))
    W[4] <- 5*sum(sp(24)) / sum(5*sp(22))
    W[6] <- 5*sum(sp(26)) / sum(5*sp(24))
    W[7] <- 5*sum(sp(27)) / sum(5*sp(25))
    W[8] <- 5*sum(sp(28)) / sum(5*sp(26))
    W[9] <- 5*sum(sp(29)) / sum(5*sp(27))
    whipple <- sum(abs(W-1), na.rm=TRUE)
  }
  return(whipple)
}


#' Whipple index (original and modified)
#' 
#' The function calculates the original and modified Whipple index to evaluate
#' age heaping.
#' 
#' The original Whipple's index is obtained by summing the number of persons in
#' the age range between 23 and 62, and calculating the ratio of reported ages
#' ending in 0 or 5 to one-fifth of the total sample. A linear decrease in the
#' number of persons of each age within the age range is assumed. Therefore,
#' low ages (0-22 years) and high ages (63 years and above) are excluded from
#' analysis since this assumption is not plausible.
#' 
#' When the digits 0 and 5 are not reported in the data, the original Whipple
#' index varies between 0 and 100, 100 if no preference for 0 or 5 is within
#' the data. When only the digits 0 and 5 are reported in the data it reaches a
#' to a maximum of 500.
#' 
#' For the modified Whipple index, age heaping is calculated for all ten digits
#' (0-9). For each digit, the degree of preference or avoidance can be
#' determined for certain ranges of ages, and the modified Whipple index then
#' is given by the absolute sum of these (indices - 1). The index is scaled between
#' 0 and 1, therefore it is 1 if all age values end with the same digit and 0 it is
#' distributed perfectly equally.
#' 
#' @name whipple
#' @param x numeric vector holding the age of persons
#' @param method \dQuote{standard} or \dQuote{modified} Whipple index.
#' @param weight numeric vector holding the weights of each person
#' @return The original or modified Whipple index.
#' @author Matthias Templ, Alexander Kowarik
#' @seealso \code{\link{sprague}}
#' @references Henry S. Shryock and Jacob S. Siegel, Methods and Materials of
#' Demography (New York: Academic Press, 1976)
#' @keywords arith
#' @export
#' @examples
#' 
#' #Equally distributed
#' age <- sample(1:100, 5000, replace=TRUE)
#' whipple(age)
#' whipple(age,method="modified")
#' 
#' # Only 5 and 10
#' age5 <- sample(seq(0,100,by=5), 5000, replace=TRUE)
#' whipple(age5)
#' whipple(age5,method="modified")
#' 
#' #Only 10
#' age10 <- sample(seq(0,100,by=10), 5000, replace=TRUE)
#' whipple(age10)
#' whipple(age10,method="modified")
#' 
whipple <- function(x, method="standard",weight=NULL){
  if(method == "standard"){
	if(is.null(weight)){
	  x <- x[x >= 23 & x <= 62]
      xm <- x %% 5
      return((length(xm[xm==0])/length(x))*500)
    }else{
	  weight <- weight[x >= 23 & x <= 62]
	  x <- x[x >= 23 & x <= 62]
	  xm <- x %% 5
	  return((sum(weight[xm==0])/sum(weight))*500)
	}
  }else if(method == "modified"){
    
	if(is.null(weight)){
		tab <- table(x)	
	}else{
		tab <- tableWt(x,weight)
	}
    W <- numeric(10)
	for(i in 1:10){
		W[i] <- sum(tab[as.numeric(names(tab))%in%seq(i-10,200,by=10)]) / (sum(tab)/10)	
	}
    return(sum(abs(W-1), na.rm=TRUE)/18)
  }else{
    stop(paste("Supplied mehtod",method,"is not implemented"))
  }
}


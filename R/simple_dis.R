#' @name simple_dis
#' @rdname simple_dis
#' @aliases univariate.dis conditional.dis
#' @title Simple generation of new variables
#' @description Fast simulation of new variables based on univariate distributions
#' @param puf data for which one additional column specified by function argument \sQuote{additional} is simulated 
#' @param data donor data
#' @param additional name of variable to be simulated
#' @param conditional conditioning variable 
#' @param weights sampling weights from data
#' @param value if \dQuote{data} then the puf including the additional variable is returned, otherwise only
#' the simulated vector.
#' @details Function uni.distribution: random draws from the weighted univariate distribution of
#' the original data ´
#' 
#' Function conditional.dis: random draws from the weighted conditional distribution
#' (conditioned on a factor variable)
#' 
#' This are simple functions to produce structural variables, variables that 
#' should have the same categories as given ones. For more advanced methods see simCategorical()
#' @seealso \code{\link{simCategorical}} 
#' @author Lydia Spies, Matthias Templ
NULL

#' @rdname simple_dis
#' @name univariate.dis
#' @examples
#' ## we don't have original data, so let's use eusilc
#' data(eusilc13)
#' data(eusilcS)
#' v1 <- univariate.dis(eusilcS, eusilc13, additional = "pl051", weights = "rb050", value = "vector")
#' table(v1)
#' table(eusilc13$pl051)
#' @export
univariate.dis <- function(puf, data, additional, weights, value = "data"){
  if(!(additional %in% colnames(data))) stop(paste("variable", additional, "not in data"))
  if(length(additional) != 1 | !is.character(additional)) stop("additional must be a string (vector of lenght 1)")
  if(length(weights) != 1 | !is.character(weights)) stop("weights must be a string (vector of lenght 1)")
  if(is.null(dim(data))) stop("data must be a matrix or data.frame")
  if(is.null(dim(puf))) stop("puf must be a matrix or data.frame")  
  if(sum(is.na(data[, additional])) > 0 & sum(is.na(data[, additional])) != dim(data)[1]) {
    var <- factorNA(data[,additional], always=TRUE)
  } else if (sum(is.na(data[, additional])) == dim(data)[1]) {
    var <- factor(c(NA, data[, additional]), exclude=c())[-1]
  } else {
    var <- as.factor(data[,additional])
  }
  tab <- tableWt(var, weights = data[, weights])
  p <- as.numeric(tab / sum(data[, weights]))
  simvar <- sample(x=levels(var)[levels(var) %in% names(tab)], size=nrow(puf), prob=p, replace=TRUE)
  if(is.factor(data[, additional])) simvar <- factor(simvar)
  if(value == "data"){
    puf[, additional] <- simvar
    return(puf[,additional]) 
  } else {
    return(simvar)
  }
}

#' @rdname simple_dis
#' @name conditional.dis
#' @examples
#' ## we don't have original data, so let's use eusilc
#' data(eusilc13)
#' data(eusilcS)
#' v1 <- conditional.dis(eusilcS, eusilc13, additional = "pl051", 
#'                            conditional = "db040", weights = "rb050")
#' table(v1) / sum(table(v1))
#' table(eusilc13$pl051) / sum(table(eusilc13$pl051))                          
#' @export
conditional.dis <- function(puf, data, additional, conditional, weights, value = "data"){
  if(!(additional %in% colnames(data))) stop(paste("variable", additional, "not in data"))
  if(!(conditional %in% colnames(data))) stop(paste("variable", conditional, "not present in data"))
  if(!(conditional %in% colnames(data))) stop(paste("variable", conditional, "not present in puf"))
  if(length(additional) != 1 | !is.character(additional)) stop("additional must be a string (vector of lenght 1)")
  if(length(conditional) != 1 | !is.character(conditional)) stop("conditional must be a string (vector of lenght 1)")
  if(length(weights) != 1 | !is.character(weights)) stop("weights must be a string (vector of lenght 1)")
  if(is.null(dim(data))) stop("data must be a matrix or data.frame")
  if(is.null(dim(puf))) stop("puf must be a matrix or data.frame")  
  if (sum(is.na(data[,additional])) > 0 & sum(is.na(data[, additional])) != dim(data)[1]) {
    var <- factorNA(data[,additional],always=TRUE)
  } else if (sum(is.na(data[,additional])) == dim(data)[1]) {
    var <- factor(c(NA, data[,additional]), exclude=c())[-1]
  } else {
    var <- as.factor(data[,additional])
  }
  ## conditional variable in data
  cond <- data[, conditional]
  ## conditional variable in puf
  condpuf <- puf[, conditional]
  ## weights variable in data
  weights <- data[, weights]
  ## levels of conditional variables in puf
  lev <- levels(puf[, conditional])
  ## levels of additional variables in data
  levadd <- levels(var)
  simvar <- numeric(nrow(puf))
  
  for (i in lev){
    group <- cond == i
    sgpuf <- sum(condpuf == i)
    tab <- tableWt(var[group], weights = weights[group])
    p <- as.numeric(tab / sum(tab))
    s <- sample(levadd[levadd %in% names(tab)],
                size = max(sum(group), sgpuf), 
                prob = p, 
                replace = TRUE)[1:sgpuf]
    simvar[condpuf == i] <- s
        
  }
  if(is.factor(data[, additional])) simvar <- factor(simvar)
  if(value == "data"){
    puf[, additional] <- simvar
    return(puf[,additional]) 
  } else {
    return(simvar)
  }
  return(puf)
}
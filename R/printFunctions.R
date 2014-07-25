setMethod("show", "dataObj", function(object){
  if ( is.null(object@ispopulation) ) {
    cat("this object does not contain any data yet!\n")
  } else {
    if(object@ispopulation){
      dname <- "population"
    } else {
      dname <- "survey sample"
    }
    cat("\n -------------- \n")
    cat(dname, "of size", nrow(object@data), "x", ncol(object@data), "\n")
    cat("\n Selected important variables: \n")
    cat("\n household ID:", object@hhid)
    cat("\n personal ID:", object@pid)
    cat("\n variable household size:", object@hhsize)
    cat("\n sampling weight:", object@weight)
    cat("\n strata:", object@strata)
    cat("\n -------------- \n")
    cat("\n")
  }
})

setMethod("show", "synthPopObj", function(object){
  if ( is.null(object@pop) ) {
    cat("synthetic population has not been generated!\n")
  } else {
    dname <- "synthetic population"
    cat("\n")
    cat("-------------- \n")
    cat(dname, " of size \n", nrow(object@pop@data), "x", ncol(object@pop@data), "\n")
    cat("\n")
    cat("build from a sample of size \n")
    cat(nrow(object@sample@data),"x",ncol(object@sample@data))
    cat("\n")
    cat("-------------- \n\n")
    cat("variables in the population:\n")
    cat(paste(colnames(object@pop@data), collapse=","))
    cat("\n")
  }
})

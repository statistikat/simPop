setMethod("show", "dataObj",
          function(object){
            if(object@ispopulation){ 
              dname <- "population"
            } else {
              dname <- "survey sample"
            }        
            cat(dname, "of size", nrow(object@data), "x", ncol(object@data), "\n")
            cat("\n household ID:", object@hhid)
            cat("\n personal ID:", object@pid)
            cat("\n variable household size:", object@hhsize)
            cat("\n sampling weight:", object@weight)
            cat("\n strata:", object@strata)
          }
)

setMethod("show", "synthPopObj",
          function(object){
              dname <- "population"
            cat(dname, " of size", nrow(object@data), "x", ncol(object@data), "\n")
            cat("\n household ID:", object@hhid)
            cat("\n personal ID:", object@pid)
            cat("\n variable household size:", object@hhsize)
            cat("\n sampling weight:", object@weight)
            cat("\n strata:", object@strata)
          }
)

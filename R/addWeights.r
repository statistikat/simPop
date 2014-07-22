setGeneric("addWeights<-", function(object, value) standardGeneric("addWeights<-"))
setReplaceMethod("addWeights", signature = c("dataObj", "list"), definition = function(object, value) {
  if ( length(value) != 2 ) {
    stop("The provided input must be the output of 'calibWeights()'\n")
  }
  if ( sum(names(value) != c("g_weights", "final_weights")) > 0 ) {
    stop("The provided input is not valid. It must be an output of 'calibWeights()'!\n")
  }
  if ( length(value[[1]]) != nrow(object@data) ) {
    stop("The dimensions do not match. Please check your input!\n")
  }
  object@data[[object@weight]] <- value$final_weights
  validObject(object)
  invisible(object)
})

setReplaceMethod("addWeights", signature = c("synthPopObj", "list"), definition = function(object, value) {
  if ( is.null(object@sample) ) {
    stop("No sample information is provided in the input object!\n")
  }
  dat <- object@sample
  addWeights(dat) <- value
  object@sample <- dat
  validObject(object)
  invisible(object)
})

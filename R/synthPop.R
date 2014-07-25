setGeneric("auxData", function(object) standardGeneric("auxData"))
setGeneric("auxData<-", function(object, value) standardGeneric("auxData<-"))

setMethod("auxData", signature="synthPopObj", definition=function(object) {
  invisible(object@sample)
})
setReplaceMethod("auxData", signature = c("synthPopObj", "dataObj"), definition = function(object, value) {
  object@sample <- value
  validObject(object)
  invisible(object)
})


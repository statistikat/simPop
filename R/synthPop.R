setGeneric("sampleObj", function(object) standardGeneric("sampleObj"))
setGeneric("sampleObj<-", function(object, value) standardGeneric("sampleObj<-"))
setMethod("sampleObj", signature="simPopObj", definition=function(object) {
  invisible(object@sample)
})
setReplaceMethod("sampleObj", signature = c("simPopObj", "dataObj"), definition = function(object, value) {
  object@sample <- value
  validObject(object)
  invisible(object)
})

setGeneric("sampleData", function(object) standardGeneric("sampleData"))
setMethod("sampleData", signature="simPopObj", definition=function(object) {
  invisible(sampleObj(object)@data)
})

setGeneric("popObj", function(object) standardGeneric("popObj"))
setGeneric("popObj<-", function(object, value) standardGeneric("popObj<-"))
setMethod("popObj", signature="simPopObj", definition=function(object) {
  invisible(object@pop)
})
setReplaceMethod("popObj", signature = c("simPopObj", "dataObj"), definition = function(object, value) {
  object@pop <- value
  validObject(object)
  invisible(object)
})

setGeneric("popData", function(object) standardGeneric("popData"))
setMethod("popData", signature="simPopObj", definition=function(object) {
  invisible(popObj(object)@data)
})

setGeneric("tableObj", function(object) standardGeneric("tableObj"))
setMethod("tableObj", signature="simPopObj", definition=function(object) {
  invisible(object@table)
})


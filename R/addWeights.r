#' Methods for function \code{addWeights}
#'
#' allows to modify sampling weights of an \code{\linkS4class{dataObj}} or
#' \code{\linkS4class{simPopObj}}-object. As input the output of
#' \code{\link{calibSample}} must be used.
#'
#' @name addWeights
#' @docType methods
#' @section Methods: \describe{ 
#' \item{list(signature(object="dataObj", value="list"))}{ \code{addWeight(inp) <- x} replaces the variable in slot
#' "data" containing sampling weights (specified by slot "weight" of input
#' \code{inp} with \code{x}. \code{x} must be a list that was computed using
#' \code{\link{calibWeights}}. } 
#' \item{list(signature(object="simPopObj", value="list"))}{ \code{addWeight(inp) <- x} updates slot "sample" of the
#' provided object of class \code{\linkS4class{simPopObj}}. It replaces the
#' variable in slot "data" that contains sampling weights (specified by slot
#' "weight") of this slot with \code{x}. \code{x} must be a list that was
#' computed using \code{\link{calibWeights}}. }}
#' @keywords methods
#' @export 
#' @examples
#' data(eusilcS)
#' data(totalsRG)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' \dontrun{
#' ## approx. 20 seconds ...
#' addWeights(inp) <- calibSample(inp, totalsRG)
#' }
setGeneric("addWeights<-", function(object, value) standardGeneric("addWeights<-"))
NULL

#' @rdname addWeights
#' @export
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


#' @rdname addWeights
#' @export
setReplaceMethod("addWeights", signature = c("simPopObj", "list"), definition = function(object, value) {
  if ( is.null(object@sample) ) {
    stop("No sample information is provided in the input object!\n")
  }
  dat <- object@sample
  addWeights(dat) <- value
  object@sample <- dat
  validObject(object)
  invisible(object)
})


#' 
#' Kish Factor
#' 
#' Compute the kish factor for a specific weight vector
#'
#' @name kishFactor
#' @param w a numeric vector with weights
#' @return The function will return the the kish factor 
#' @author Alexander Kowarik
#' @export kishFactor
#' @examples
#' kishFactor(rep(1,10))
#' kishFactor(rlnorm(10))
kishFactor <- function(w){
  
  if(!requireNamespace("surveysd",quietly=TRUE)){
    stop("the functionality of kishFactor moved to survey sd and the package survey sd is not available")
  }
  .Deprecated("kishFactor", package="surveysd",
              msg="the function kishFactor moved to the package surveysd in the function kishFactor")
  return(surveysd::kishFactor(w))
}
#' Iterative Proportional Updating
#' 
#' Adjust sampling weights to given totals based on household-level and/or
#' individual level constraints.
#' 
#' This function implements the weighting procedure described 
#' [here](http://www.ajs.or.at/index.php/ajs/article/viewFile/doi10.17713ajs.v45i3.120/512). 
#' 
#' `conP` and `conH` are contingency tables, which can be created with `xtabs`. The `dimnames` of those
#' tables should match the names and levels of the corresponding columns in `dat`.
#' 
#' `maxIter`, `epsP` and `epsH` are the stopping criteria. `epsP` and `epsH` describe relative tolerances
#' in the sense that
#' \out{\deqn{1-epsP < \frac{w_{i+1}}{w_i} < 1+epsP}{1-epsP < w(i+1)/w(i) < 1+epsP} }
#' will be used as convergence criterium. Here i is the iteration step and wi is the weight of a 
#' specific person at step i.
#' 
#' The algorithm 
#' performs best if all varables occuring in the constraints (\code{conP} and \code{conH}) as well as the 
#' household variable are coded as \code{factor}-columns in \code{dat}. Otherwise, conversions will be necessary
#' which can be monitored with the `conversion_messages` argument.
#' Setting `check_hh_vars` to `FALSE` can also incease the performance of the scheme.
#'
#' @name ipu2 
#' @md
#' @aliases ipu2
#' @param dat a \code{data.table} containing household ids (optionally), base
#' weights (optionally), household and/or personal level variables (numerical
#' or categorical) that should be fitted.
#' @param hid name of the column containing the household-ids
#' within \code{dat} or NULL if such a variable does not exist.
#' @param w name if the column containing the base
#' weights within \code{dat} or NULL if such a variable does not exist. In the
#' latter case, every observation in \code{dat} is assigned a starting weight
#' of 1.
#' @param conP list or (partly) named list defining the constraints on person
#' level.  The list elements are contingency tables in array representation
#' with dimnames corresponding to the names of the relevant calibration
#' variables in \code{dat}. If a numerical variable is to be calibrated, the
#' respective list element has to be named with the name of that numerical
#' variable. Otherwise the list element shoud NOT be named.
#' @param conH list or (partly) named list defining the constraints on
#' household level.  The list elements are contingency tables in array
#' representation with dimnames corresponding to the names of the relevant
#' calibration variables in \code{dat}. If a numerical variable is to be
#' calibrated, the respective list element has to be named with the name of
#' that numerical variable. Otherwise the list element shoud NOT be named.
#' @param epsP numeric value or list (of numeric values and/or arrays)
#' specifying the convergence limit(s) for \code{conP}. The list can contain
#' numeric values and/or arrays which must appear in the same order as the
#' corresponding constraints in \code{conP}. Also, an array must have the same
#' dimensions and dimnames as the corresponding constraint in \code{conP}.
#' @param epsH numeric value or list (of numeric values and/or arrays)
#' specifying the convergence limit(s) for \code{conH}. The list can contain
#' numeric values and/or arrays which must appear in the same order as the
#' corresponding constraints in \code{conH}. Also, an array must have the same
#' dimensions and dimnames as the corresponding constraint in \code{conH}.
#' @param verbose if TRUE, some progress information will be printed.
#' @param bound numeric value specifying the multiplier for determining the
#' weight trimming boundary if the change of the base weights should be
#' restricted, i.e. if the weights should stay between 1/\code{bound}*\code{w}
#' and \code{bound}*\code{w}.
#' @param maxIter numeric value specifying the maximum number of iterations
#' that should be performed.
#' @param meanHH if TRUE, every person in a household is assigned the mean of
#' the person weights corresponding to the household. If \code{"geometric"}, the geometric mean
#' is used rather than the arithmetic mean.
#' @param allPthenH if TRUE, all the person level calibration steps are performed before the houshold level calibration steps (and \code{meanHH}, if specified). 
#' If FALSE, the houshold level calibration steps (and \code{meanHH}, if specified) are performed after everey person level calibration step.
#' This can lead to better convergence properties in certain cases but also means that the total number of calibration steps is increased.
#' @param returnNA if TRUE, the calibrated weight will be set to NA in case of no convergence.
#' @param looseH if FALSE, the actual constraints \code{conH} are used for calibrating all the hh weights. 
#' If TRUE, only the weights for which the lower and upper thresholds defined by \code{conH} and \code{epsH} are exceeded
#' are calibrated. They are however not calibrated against the actual constraints \code{conH} but against
#' these lower and upper thresholds, i.e. \code{conH}-\code{conH}*\code{epsH} and \code{conH}+\code{conH}*\code{epsH}.
#' @param numericalWeighting If NULL computeLinear from the pacakge survey sd will be used.
#' @param check_hh_vars If \code{TRUE} check for non-unique values inside of a household for variables in 
#'                      household constraints
#' @param conversion_messages show a message, if inputs need to be reformatted. This can be useful for speed 
#'        optimizations if ipu2 is called several times with similar inputs (for example bootstrapping)
#' @return The function will return the input data \code{dat} with the
#' calibrated weights \code{calibWeight} as an additional column as well as attributes. If no convergence has been reached in `maxIter` 
#' steps, and `returnNA` is `TRUE` (the default), the column `calibWeights` will only consist of `NA`s. The attributes of the table are
#' attributes derived from the `data.table` class as well as the following.
#' \tabular{ll}{
#'   `converged` \tab Did the algorithm converge in `maxIter` steps? \cr
#'   `iterations` \tab The number of iterations performed. \cr
#'   `conP`, `conH`, `epsP`, `epsH` \tab See Arguments. \cr
#'   `conP_adj`, `conH_adj` \tab Adjusted versions of `conP` and `conH` \cr
#'   `formP`, `formH` \tab Formulas that were used to calculate `conP_adj` and `conH_adj` based on the output table.
#' }
#' @seealso \code{\link{ipu}}
#' @export ipu2
#' @author Alexander Kowarik, Gregor de Cillia
#' @examples
#' library(data.table)
#' data(eusilcS)
#' setDT(eusilcS)
#' eusilcS <- eusilcS[, list(db030,hsize,db040,age,rb090,netIncome,db090,rb050)]
#' 
#' ## rename columns
#' setnames(eusilcS, "rb090", "gender")
#' setnames(eusilcS, "db040", "state")
#' setnames(eusilcS, "db030", "household")
#' setnames(eusilcS, "rb050", "weight")
#' 
#' ## some recoding
#' # generate age groups
#' eusilcS[, agegroup := cut(age, c(-Inf, 10*1:9, Inf), right = FALSE)]
#' # some recoding of netIncome for reasons of simplicity
#' eusilcS[is.na(netIncome), netIncome := 0] 
#' eusilcS[netIncome < 0, netIncome := 0] 
#' # set hsize to 1,...,5+
#' eusilcS[, hsize := cut(hsize, c(0:4, Inf), labels = c(1:4, "5+"))]
#' # treat households as a factor variable
#' eusilcS[, household := as.factor(household)]
#' 
#' ## example for base weights assuming a simple random sample of households stratified per region
#' eusilcS[, regSamp := .N, by = state]
#' eusilcS[, regPop := sum(weight), by = state]
#' eusilcS[, baseWeight := regPop/regSamp]
#' 
#' ## constraints on person level
#' # age 
#' conP1 <- xtabs(weight ~ agegroup, data = eusilcS)
#' # gender by region
#' conP2 <- xtabs(weight ~ gender + state, data = eusilcS)
#' # personal net income by gender
#' conP3 <- xtabs(weight*netIncome ~ gender, data = eusilcS)
#' 
#' ## constraints on household level
#' conH1 <- xtabs(weight ~ hsize + state, data = eusilcS, subset = !duplicated(household))
#' 
#' # array of convergence limits for conH1
#' epsH1 <- conH1
#' epsH1[1:4,] <- 0.005
#' epsH1["5+",] <- 0.2
#' 
#' # without array epsH1
#' calibweights1 <- ipu2(eusilcS, hid = "household", 
#'                       conP = list(conP1, conP2, netIncome = conP3), 
#'                       conH = list(conH1), 
#'                       epsP = list(1e-06, 1e-06, 1e-03),
#'                       epsH = 0.01,  
#'                       bound = NULL, verbose = TRUE,  maxIter = 200)
#' 
#' # with array epsH1, base weights and bound
#' calibweights2 <- ipu2(eusilcS, hid = "household", 
#'                       conP = list(conP1, conP2), 
#'                       conH = list(conH1), 
#'                       epsP = 1e-06,
#'                       epsH = list(epsH1),  
#'                       w = "baseWeight",
#'                       bound = 4, verbose = TRUE, maxIter = 200)
#'                       
#' # show an adjusted version of conP and the original
#' attr(calibweights2, "conP_adj")
#' attr(calibweights2, "conP")
#' 
ipu2 <- function(dat,hid=NULL,conP=NULL,conH=NULL,epsP=1e-6,epsH=1e-2,verbose=FALSE,
                 w=NULL,bound=4,maxIter=200,meanHH=TRUE,allPthenH=TRUE,returnNA=TRUE,looseH=FALSE,
                 numericalWeighting=NULL, check_hh_vars = TRUE, conversion_messages = FALSE){
  .Deprecated("ipf",package="surveysd", msg="the functionality of ipu2 moved to the package surveysd in the function ipf")
  if(!requireNamespace("surveysd",quietly=TRUE)){
    stop("the functionality of ipu2 moved to surveysd::ipf and the package surveysd is not available")
  }
  if(is.null(numericalWeighting)){
    numericalWeighting <- surveysd::computeLinear
  }
  return(surveysd::ipf(dat, hid, conP, conH, epsP, epsH, verbose, w, bound, maxIter, 
              meanHH, allPthenH, returnNA, looseH, numericalWeighting,
              check_hh_vars, conversion_messages))
  
}

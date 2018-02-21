getMeanFun <- function(meanHH){
  if(isTRUE(meanHH))
    meanHH <- "arithmetic"
  if(identical(meanHH, FALSE))
    meanHH <- "none"
  meanfun <- switch (meanHH,
    arithmetic = arithmetic_mean,
    geometric = geometric_mean,
    none = function(x,w){x}
  )
  if(is.null(meanfun))
    stop("invalid value for meanHH")
  meanfun
}

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
  if(!is.numeric(w)){
    stop("The input must be a numeric vector")
  }
  n <- length(w)
  sqrt(n*sum(w^2)/sum(w)^2)  
}
boundsFak <- function(g1,g0,f,bound=4){ # Berechnet die neuen Gewichte (innerhalb 4, .25 Veraenderungsraten)
  g1 <- g1 * f
  TF <- which((g1/g0)>bound)
  TF[is.na(TF)] <- FALSE
  g1[TF] <- bound*g0[TF]
  TF <- which((g1/g0)<(1/bound))
  TF[is.na(TF)] <- FALSE
  g1[TF] <- (1/bound)*g0[TF]
  return(g1)
}
boundsFakHH <- function(g1,g0,eps,orig,p,bound=4){ # Berechnet die neuen Gewichte fuer Unter- und Obergrenze (innerhalb 4, .25 Veraenderungsraten)
  u <- orig*(1 - eps)
  o <- orig*(1 + eps)
  
  pbo <- which(p > o)
  psu <- which(p < u)
  g1[pbo] <- g1[pbo] * o[pbo]/p[pbo]
  g1[psu] <- g1[psu] * u[psu]/p[psu]

  TF <- which((g1/g0)>bound)
  g1[TF] <- bound*g0[TF]
  TF <- which((g1/g0)<(1/bound))
  g1[TF] <- (1/bound)*g0[TF]
  return(g1)
}

check_population_totals <- function(con){
  # do not apply this check for numerical calibration
  ind <- which(names(con) == "")
  
  if (length(ind) == 0)
    return(NULL)
  
  pop_totals <- vapply(ind, function(index){ sum(con[[index]]) }, 0)
  rel_errors <- abs(pop_totals - pop_totals[1])/pop_totals[1]
  
  # use a 1% tolerance. Maybe it would be better to make the tolerance dependent on conH?
  if(any(rel_errors > 1e-2))
    stop("population totals for different constraints do not match")
}

calibP <- function(i,conP, epsP, dat, error, valueP, pColNames, bound, verbose, calIter, numericalWeighting){
  OriginalSortingVariable <- V1 <- baseWeight <- calibWeight <- epsvalue <- f <- NULL
  temporary_hid <- temporary_hvar <- tmpVarForMultiplication <- value <- wValue <- wvst<- NULL
  
  if(is.list(epsP)){
    epsPcur <- epsP[[i]]
  }else{
    epsPcur <- epsP
  }
  combined_factors <- dat[[paste0("combined_factors_", i)]]
  
  if(isTRUE(names(conP)[i]!="")){
    ## numerical variable to be calibrated
    ## use name of conP list element to define numerical variable
    setnames(dat,names(conP)[i],"tmpVarForMultiplication")
    setnames(dat,valueP[i],"value")
    
    dat[, f := ipu_step_f(calibWeight*tmpVarForMultiplication, 
                          combined_factors, conP[[i]])]
    dat[, wValue := value/f]
    
    # try to divide the weight between units with larger/smaller value in the numerical variable linear
    dat[,f:=numericalWeighting(head(wValue,1),head(value,1),tmpVarForMultiplication,calibWeight),
        by=combined_factors]
    
    if(is.array(epsPcur)){## for numeric variables not the factor f is used but the abs relative deviation is computed per class
      curEps <- abs(dat[!is.na(f),1/(value/wValue)-1]) ## curEps is computed for all observations to compare it with the right epsValue
    }else{
      curEps <- dat[!is.na(f),max(abs(1/(value/wValue)-1))]  ## only the max curEps is needed because it is the same for all classed of the current variable
    }
    
    setnames(dat,"tmpVarForMultiplication",names(conP)[i])
    setnames(dat,"value",valueP[i])
  }else{
    # categorical variable to be calibrated
    dat[, f := ipu_step_f(dat$calibWeight, combined_factors, conP[[i]])]
    
    if(is.array(epsPcur)){
      curEps <- abs(dat[!is.na(f),1/f-1]) ## curEps is computed for all observations to compare it with the right epsValue
    }else{
      curEps <- dat[!is.na(f),max(abs(1/f-1))]  ## only the max curEps is needed because it is the same for all classed of the current variable
    }
  }
  if(is.array(epsPcur)){
    epsPcur <- dat[!is.na(f), paste0("epsP_", i)]
  }
  
  if(any(curEps>epsPcur)){## sicherheitshalber abs(epsPcur)? Aber es wird schon niemand negative eps Werte uebergeben??
    
    if(!is.null(bound)){
      dat[!is.na(f),calibWeight:=boundsFak(calibWeight,baseWeight,f,bound=bound)]#,by=eval(pColNames[[i]])]  
    }else{
      dat[!is.na(f),calibWeight:=f*calibWeight,by=combined_factors]
    }
    error <- TRUE
  }
  
  if(verbose&&any(curEps>epsPcur)&&calIter%%10==0){
    if(calIter%%100==0)
      print(subset(dat,!is.na(f))[abs(1/f-1)>epsPcur][,list(mean(f)-1,.N),by=eval(pColNames[[i]])])
    #print(dat[abs(1/f-1)>epsPcur][,list(mean(f),.N),by=eval(pColNames[[i]])])
    cat(calIter, ":Not yet converged for P-Constraint",i,"\n")
  }
  
  return(error)
}

calibH <- function(i,conH, epsH, dat, error, valueH, hColNames, bound, verbose, calIter, looseH){  
  OriginalSortingVariable <- V1 <- baseWeight <- calibWeight <- epsvalue <- f <- NULL
  temporary_hid <- temporary_hvar <- tmpVarForMultiplication <- value <- wValue <- wvst<- NULL
  
  if(is.list(epsH)){
    epsHcur <- epsH[[i]]
  }else{
    epsHcur <- epsH
  }
  
  setnames(dat,valueH[i],"value")
  combined_factors <- dat[[paste0("combined_factors_h_", i)]]
  
  if(isTRUE(names(conH)[i]!="")){
    ## numerical variable to be calibrated
    ## use name of conH list element to define numerical variable
    setnames(dat,names(conH)[i],"tmpVarForMultiplication")
    
    dat[, f := ipu_step_f(calibWeight*wvst*tmpVarForMultiplication, 
                          combined_factors, conH[[i]])]
    
    setnames(dat,"tmpVarForMultiplication",names(conH)[i])
  }else{
    # categorical variable to be calibrated
    dat[, f := ipu_step_f(calibWeight*wvst, combined_factors, conH[[i]])]
  }
  
  dat[, wValue := value/f]
  
  if(is.array(epsHcur)){
    curEps <- dat[,abs(1/f-1)]
    epsHcur <- dat[[paste0("epsH_", i)]]
  }else{
    curEps <- dat[,max(abs(1/f-1))]
  }

  if(any(curEps>epsHcur)){    
    if(!is.null(bound)){
      if(!looseH){
        dat[,calibWeight:=boundsFak(g1=calibWeight,g0=baseWeight,f=f,bound=bound)]#,by=eval(hColNames[[i]])]    
      }else{
        dat[,calibWeight:=boundsFakHH(g1=calibWeight,g0=baseWeight,eps=epsHcur*.9,orig=value,p=wValue,bound=bound)]  
      }
    }else{
      dat[,calibWeight:=f*calibWeight,by=combined_factors]
    }
    error <- TRUE
  }
  
  setnames(dat,"value",valueH[i])
  
  if(verbose&&any(curEps>epsHcur)&&calIter%%10==0){
    if(calIter%%100==0)
      print(subset(dat,!is.na(f))[abs(1/f-1)>epsHcur][,list(mean(f)-1,.N),by=eval(hColNames[[i]])])
    cat(calIter, ":Not yet converged for H-Constraint",i,"\n")
  }
  return(error)
}

## recreate the formula argument to xtabs based on conP, conH
getFormulas <- function(con, w = "calibWeight"){
  formOut <- NULL
  for(i in seq_along(con)){
    lhs <- names(con)[i]
    if(is.null(lhs) || lhs == ""){
      lhs = w
    }
    rhs <- paste(names(dimnames(con[[i]])), collapse = "+")
    formOut[[i]] <- formula(paste(lhs, "~", rhs), env = .GlobalEnv)
  }
  formOut
}

## enrich dat_original with the calibrated weights and assign attributes
addWeightsAndAttributes <- function(dat, conP, conH, epsP, epsH, dat_original, maxIter, calIter, returnNA){
  outTable <- copy(dat_original)
  
  # add calibrated weights. Use setkey to make sure the indexes match
  setkey(dat, OriginalSortingVariable)
  outTable[, calibWeight := ifelse((maxIter < calIter) & returnNA, NA, dat$calibWeight)]
  
  formP <- getFormulas(conP)
  formH <- getFormulas(conH)

  # general information  
  setattr(outTable, "converged", (maxIter >= calIter))
  setattr(outTable, "iterations", min(maxIter, calIter)) # return maxIter in case of no convergence

  # input constraints
  setattr(outTable, "conP", conP)
  setattr(outTable, "conH", conH)
  
  # adjusted constraints (conP, conH according to the calibrated weights)
  setattr(outTable, "conP_adj", lapply(formP, xtabs, dat))
  setattr(outTable, "conH_adj", lapply(formH, xtabs, dat))
  
  # tolerances
  setattr(outTable, "epsP", epsP)
  setattr(outTable, "epsH", epsH)
  
  # formulas
  setattr(outTable, "formP", formP)
  setattr(outTable, "formH", formH)
  
  # not used yet
  #class(outTable) <- c("ipu2", class(outTable))
  
  invisible(outTable)
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
#' @param numericalWeighting See [numericalWeighting]
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
                 numericalWeighting=computeLinear, check_hh_vars = TRUE, conversion_messages = FALSE){
  
  check_population_totals(conP)
  check_population_totals(conH)

  OriginalSortingVariable <- V1 <- baseWeight <- calibWeight <- epsvalue <- f <- NULL
  temporary_hid <- temporary_hvar <- tmpVarForMultiplication <- value <- wValue <- wvst<- NULL
  dat_original <- dat
  dat <- copy(dat)
  ## originalsorting is fucked up without this
  dat[,OriginalSortingVariable:=.I]
  meanfun <- getMeanFun(meanHH)
  
  # dat sollte ein data.table sein
  # w ein Name eines Basisgewichts oder NULL
  valueP <- paste0("valueP",seq_along(conP))###fixed target value, should not be changed in iterations
  valueH <- paste0("valueH",seq_along(conH))
  ###Housekeeping of the varNames used
  usedVarNames <- c(valueP,valueH,"value","baseWeight","wvst","wValue")

  if(any(names(dat)%in%usedVarNames)){
    renameVars <- names(dat)[names(dat)%in%usedVarNames]
    setnames(dat,renameVars,paste0(renameVars,"_safekeeping"))
    if(isTRUE(w=="baseWeight"))
      w <- "baseWeight_safekeeping"
  }
  ### Treatment of HID, creating 0,1 var for being the first hh member
  delVars <- c()
  if(is.null(hid)){
    delVars <- c("hid")
    hid <- "hid"
    dat[,hid:=1:nrow(dat)]
    dat[,wvst:=1] 
  }else{
    setnames(dat,hid,"temporary_hid")
    dat[,wvst:=as.numeric(!duplicated(temporary_hid))]
    setnames(dat,"temporary_hid",hid)
  }
  
  setnames(dat, hid, "temporary_hid")
  if(!is.factor(dat$temporary_hid)){
    if(conversion_messages)
      message("convert household variable ", hid, " to factor")
    dat[, temporary_hid := as.factor(temporary_hid)]
  }
  setnames(dat, "temporary_hid", hid)
  
  ## Names of the calibration variables for Person and household dimension
  pColNames <- lapply(conP,function(x)names(dimnames(x)))
  hColNames <- lapply(conH,function(x)names(dimnames(x)))
  
  for(i in seq_along(conP)){
    current_colnames <- pColNames[[i]]

    for(colname in current_colnames){
      if(!inherits(dat[[colname]], "factor")){
        if(conversion_messages)
          message("converting column ", colname, " to factor")
        set(
          dat, j = colname, 
          value = factor(dat[[colname]], levels = dimnames(conP[[i]])[[colname]])
        )
      }
      else if(!identical(levels(dat[[colname]]), dimnames(conP[[i]])[[colname]])){
        if(conversion_messages)
          message("correct levels of column ", colname)
        set(
          dat, j = colname, 
          value = factor(dat[[colname]], levels = dimnames(conP[[i]])[[colname]])
        )
      }
    }
    combined_factors <- combine_factors(dat, pColNames[[i]])
    
    dat[, paste0("combined_factors_", i) := combined_factors]
    dat[, paste0("valueP", i) := conP[[i]][combined_factors]]
  }
  for(i in seq_along(conH)){
    colnames <- hColNames[[i]]
    
    ## make sure the columns mentioned in the contingency table are in fact factors
    for(colname in colnames){
      if (!inherits(dat[[colname]], "factor")){
        if(conversion_messages)
          message("converting column ", colname, " to factor")
        set(
          dat, j = colname, 
          value = factor(dat[[colname]], levels = dimnames(conH[[i]])[[colname]])
        )
      }
      else if(!identical(levels(dat[[colname]]), dimnames(conH[[i]])[[colname]])){
        if(conversion_messages)
          message("correct levels of column ", colname)
        set(
          dat, j = colname, 
          value = factor(dat[[colname]], levels = dimnames(conH[[i]])[[colname]])
        )
      }
    }
    
    combined_factors <- combine_factors(dat, hColNames[[i]])
    
    dat[, paste0("combined_factors_h_", i) := combined_factors]
    dat[, paste0("valueH", i) := conH[[i]][combined_factors]]
  }
  
  if(is.null(w)){
    if(!is.null(bound)&&is.null(w))
      stop("Bounds are only reasonable if base weights are provided")
    dat[,calibWeight:=1]
    #delVars <- c(delVars,"baseWeight")
  }else{
    dat[,calibWeight:=dat[,w,with=FALSE]]
    setnames(dat,w,"baseWeight")
  }
  
  if(check_hh_vars){
    ## Check for non-unqiue values inside of a household for variabels used in Household constraints
    for(hh in hColNames){
      setnames(dat,hid,"temporary_hid")
      for(h in hh){
        setnames(dat,h,"temporary_hvar")
        if(dat[,length(unique(temporary_hvar)),by=temporary_hid][,any(V1!=1)]){
          stop(paste(h,"has different values inside a household"))
        }
        setnames(dat,"temporary_hvar",h)
      }
      setnames(dat,"temporary_hid",hid)
    }
  }
  
  if(is.list(epsP)){
    for(i in seq_along(epsP)){
      if(is.array(epsP[[i]])){
        combined_factors <- dat[[paste0("combined_factors_", i)]]
        dat[, paste0("epsP_", i) := epsP[[i]][combined_factors] ]
      }
    }
  }
  if(is.list(epsH)){
    for(i in seq_along(epsH)){
      if(is.array(epsH[[i]])){
        combined_factors <- dat[[paste0("combined_factors_h_", i)]]
        dat[, paste0("epsH_", i) := epsH[[i]][combined_factors] ]
      }
    }
  }
  
  ###Calib
  error <- TRUE
  calIter <- 1
  
  while(error&&calIter<=maxIter){
    error <- FALSE
    
    if(allPthenH){
      ### Person calib
      for(i in seq_along(conP)){
        
        error <- calibP(i=i, conP=conP, epsP=epsP, dat=dat, error=error,
                      valueP=valueP, pColNames=pColNames,bound=bound, verbose=verbose, calIter=calIter, numericalWeighting=numericalWeighting)
      }
      
      ## replace person weight with household average
      dh <- dat[[hid]]
      dat[,calibWeight := meanfun(calibWeight, dh)]

      ### Household calib
      for(i in seq_along(conH)){
        error <- calibH(i=i, conH=conH, epsH=epsH, dat=dat, error=error,
                      valueH=valueH, hColNames=hColNames,bound=bound, verbose=verbose, calIter=calIter, looseH=looseH)
      }
    }else{
      ### Person calib
      for(i in seq_along(conP)){
        
        error <- calibP(i=i, conP=conP, epsP=epsP, dat=dat, error=error,
                      valueP=valueP, pColNames=pColNames,bound=bound, verbose=verbose, calIter=calIter, numericalWeighting=numericalWeighting)

        ## replace person weight with household average
        dh <- dat[[hid]]
        dat[,calibWeight := meanfun(calibWeight, dh)]

        ### Household calib
        for(i in seq_along(conH)){
          error <- calibH(i=i, conH=conH, epsH=epsH, dat=dat, error=error,
                        valueH=valueH, hColNames=hColNames,bound=bound, verbose=verbose, calIter=calIter, looseH=looseH)
        }
      }
    }
    
    if(verbose&&!error){
      cat("Convergence reached in ",calIter," steps \n")
    }else if(verbose&&maxIter==calIter){
      cat("Not converged in",maxIter,"steps \n")
    }
    calIter <- calIter + 1 
  }
  
  addWeightsAndAttributes(dat, conP, conH, epsP, epsH, dat_original, maxIter, calIter, returnNA)
}

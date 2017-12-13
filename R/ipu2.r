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
  TF <- (g1/g0)>bound
  TF[is.na(TF)] <- FALSE
  g1[TF] <- bound*g0[TF]
  TF <- (g1/g0)<(1/bound)
  TF[is.na(TF)] <- FALSE
  g1[TF] <- (1/bound)*g0[TF]
  return(g1)
}
boundsFakHH <- function(g1,g0,eps,orig,p,bound=4){ # Berechnet die neuen Gewichte fuer Unter- und Obergrenze (innerhalb 4, .25 Veraenderungsraten)
  u <- orig*(1-eps)
  o <- orig*(1+eps)
  g1[p>o] <- g1[p>o] * o[p>o]/p[p>o]
  g1[p<u] <- g1[p<u] * u[p<u]/p[p<u]
  g1[(g1/g0)>bound] <- bound*g0[(g1/g0)>bound]
  g1[(g1/g0)<(1/bound)] <- (1/bound)*g0[(g1/g0)<(1/bound)]
  return(g1)
}
meltepsfun <- function(x){
  if(is.array(x)){
    x <- melt(x,as.is=TRUE,value.name="epsvalue")
  }
  return(x)
}
calibP <- function(i,conP, epsP, dat, error, valueP, pColNames, bound, verbose, calIter, numericalWeighting){  
  if(is.list(epsP)){
    epsPcur <- epsP[[i]]
  }else{
    epsPcur <- epsP
  }
  
  if(isTRUE(names(conP)[i]!="")){
    ## numerical variable to be calibrated
    ## use name of conP list element to define numerical variable
    setnames(dat,names(conP)[i],"tmpVarForMultiplication")
    combined_factors <- dat[[paste0("combined_factors_", i)]]
    dat[, f := ipu_step_f(calibWeight*tmpVarForMultiplication, 
                          combined_factors, conP[[i]])]
    setnames(dat,valueP[i],"value")
    dat[, wValue := f*value]
    
    # try to divide the weight between units with larger/smaller value in the numerical variable linear
    dat[,f:=numericalWeighting(head(wValue,1),head(value,1),tmpVarForMultiplication,calibWeight),by=eval(pColNames[[i]])]
    
    # if(meanHH){
    #   # Apply person-level adjustments in a multilicative way per http://www.scag.ca.gov/Documents/PopulationSynthesizerPaper_TRB.pdf
    #   # dat[duplicated(subset(dat, select = c(hid, eval(pColNames[[i]])))),f:= 1]
    #   # dat[,f:=prod(f),by=eval(hid)] 
    #   dat[,f:=gm_mean(f),by=eval(hid)] 
    # }
    
    setnames(dat,"tmpVarForMultiplication",names(conP)[i])
    if(is.array(epsPcur)){## for numeric variables not the factor f is used but the abs relative deviation is computed per class
      curEps <- abs(dat[!is.na(f),1/(value/wValue)-1]) ## curEps is computed for all observations to compare it with the right epsValue
    }else{
      curEps <- dat[!is.na(f),max(abs(1/(value/wValue)-1))]  ## only the max curEps is needed because it is the same for all classed of the current variable
    }
  }else{
    # categorical variable to be calibrated
    setnames(dat,valueP[i],"value")
    combined_factors <- dat[[paste0("combined_factors_", i)]]
    dat[, f := ipu_step_f(dat$calibWeight, combined_factors, conP[[i]])]
    
    # if(meanHH){
    #   # Apply person-level adjustments in a multilicative way per http://www.scag.ca.gov/Documents/PopulationSynthesizerPaper_TRB.pdf
    #   # dat[duplicated(subset(dat, select = c(hid, eval(pColNames[[i]])))),f:= 1]
    #   # dat[,f:=prod(f),by=eval(hid)] 
    #   dat[,f:=gm_mean(f),by=eval(hid)] 
    # }
    
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
      dat[!is.na(f),calibWeight:=f*calibWeight,by=eval(pColNames[[i]])]
    }
    error <- TRUE
  }
  setnames(dat,"value",valueP[i])
  
  if(verbose&&any(curEps>epsPcur)&&calIter%%10==0){
    if(calIter%%100==0)
      print(subset(dat,!is.na(f))[abs(1/f-1)>epsPcur][,list(mean(f),.N),by=eval(pColNames[[i]])])
    #print(dat[abs(1/f-1)>epsPcur][,list(mean(f),.N),by=eval(pColNames[[i]])])
    cat(calIter, ":Not yet converged for P-Constraint",i,"\n")
  }
  
  return(list(dat=dat,error=error))
}
calibH <- function(i,conH, epsH, dat, error, valueH, hColNames, bound, verbose, calIter, looseH){  
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
    
    dat[, f := ipu_step_f(dat$calibWeight*dat$wvst*dat$tmpVarForMultiplication, 
                          combined_factors, conH[[i]])]
    dat[, wValue := f*value]
    
    setnames(dat,"tmpVarForMultiplication",names(conH)[i])
  }else{
    # categorical variable to be calibrated
    dat[, f := ipu_step_f(dat$calibWeight*dat$wvst, combined_factors, conH[[i]])]
    dat[, wValue := f*value]
  }
  
  if(is.array(epsHcur)){
    curEps <- dat[,abs(1/f-1)]
    epsHcur <- dat[[paste0("epsH_", i)]]
  }else{
    curEps <- dat[,max(abs(1/f-1))]
  }
  # if(is.array(epsHcur)){## was soll das?
  #   dat[,epsvalue:=NULL]
  # }
  if(any(curEps>epsHcur)){    
    if(!is.null(bound)){
      if(!looseH){
        dat[,calibWeight:=boundsFak(g1=calibWeight,g0=baseWeight,f=f,bound=bound)]#,by=eval(hColNames[[i]])]    
      }else{
        dat[,calibWeight:=boundsFakHH(g1=calibWeight,g0=baseWeight,eps=epsHcur,orig=value,p=wValue,bound=bound)]  
      }
    }else{
      dat[,calibWeight:=f*calibWeight,by=eval(hColNames[[i]])]
    }
    error <- TRUE
  }
  if("epsvalue"%in%colnames(dat)){
    dat[,epsvalue:=NULL]
  }
  
  setnames(dat,"value",valueH[i])
  # if(any(curEps>epsHcur)){ 
  #   error <- TRUE
  # }
  if(verbose&&any(curEps>epsHcur)&&calIter%%10==0){
    if(calIter%%100==0)
      print(subset(dat,!is.na(f))[abs(1/f-1)>epsHcur][,list(mean(f),.N),by=eval(hColNames[[i]])])
    cat(calIter, ":Not yet converged for H-Constraint",i,"\n")
  }
  return(list(dat=dat,error=error))
}
# From package robCompositions
# gm_mean <- function(x){
#   if (!is.numeric(x)) 
#     stop("x has to be a vector of class numeric")
#   if (any(na.omit(x == 0))) 
#     0
#   else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
# }

#' Iterative Proportional Updating
#' 
#' Adjust sampling weights to given totals based on household-level and/or
#' individual level constraints.
#' 
#'
#' @name ipu2 
#' @aliases ipu2 computeLinear computeFrac 
#' @param dat a \code{data.table} containing household ids (optionally), base
#' weights (optionally), household and/or personal level variables (numerical
#' or categorical) that should be fitted.
#' @param hid character vector specifying the variable containing household-ids
#' within \code{dat} or NULL if such a variable does not exist.
#' @param w character vector specifying the variable containing the base
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
#' the person weights corresponding to the household.
#' @param allPthenH if TRUE, all the person level calibration steps are performed before the houshold level calibration steps (and \code{meanHH}, if specified). 
#' If FALSE, the houshold level calibration steps (and \code{meanHH}, if specified) are performed after everey person level calibration step.
#' This can lead to better convergence properties in certain cases but also means that the total number of calibration steps is increased.
#' @param returnNA if TRUE, the calibrated weight will be set to NA in case of no convergence.
#' @param looseH if FALSE, the actual constraints \code{conH} are used for calibrating all the hh weights. 
#' If TRUE, only the weights for which the lower and upper thresholds defined by \code{conH} and \code{epsH} are exceeded
#' are calibrated. They are however not calibrated against the actual constraints \code{conH} but against
#' these lower and upper thresholds, i.e. \code{conH}-\code{conH}*\code{epsH} and \code{conH}+\code{conH}*\code{epsH}.
#' @param numericalWeighting ...
#' @param curValue the current value of the group total
#' @param Value the target group total
#' @param numericVar vector with the values of the numeric variable
#' @param weightVec vector with the current weights
#' @param boundLinear the result of computeLinear will be bound by \code{1/boundLinear} and \code{boundLinear}
#' @param check_hh_vars If \code{TRUE} check for non-unique values inside of a household for variables in 
#'                      household constraints
#' @return The function will return the input data \code{dat} with the
#' calibrated weights \code{calibWeight} as an additional column.
#' @seealso \code{\link{ipu}}
#' @export ipu2
#' @author Alexander Kowarik
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
#' # without array epsP1
#' calibweights1 <- ipu2(eusilcS, hid = "household", 
#'                       conP = list(conP1, conP2, netIncome = conP3), 
#'                       conH = list(conH1), 
#'                       epsP = list(1e-06, 1e-06, 1e-03),
#'                       epsH = 0.01,  
#'                       bound = NULL, verbose = TRUE,  maxIter = 200)
#' 
#' # with array epsP1, base weights and bound
#' calibweights2 <- ipu2(eusilcS, hid = "household", 
#'                       conP = list(conP1, conP2), 
#'                       conH = list(conH1), 
#'                       epsP = 1e-06,
#'                       epsH = list(epsH1),  
#'                       w = "baseWeight",
#'                       bound = 4, verbose = TRUE, maxIter = 200)
#   fn <- function(a){
#    f <- a[1]*var+a[2]
#    (sum(f*var*w)-v)^2+(sum(f*w)-sum(w))^2
#  }
#  coef <- optim(c(1,1),fn)$par
ipu2 <- function(dat,hid=NULL,conP=NULL,conH=NULL,epsP=1e-6,epsH=1e-2,verbose=FALSE,
                 w=NULL,bound=4,maxIter=200,meanHH=TRUE,allPthenH=TRUE,returnNA=TRUE,looseH=FALSE,
                 numericalWeighting=computeLinear, check_hh_vars = TRUE){
  
  OriginalSortingVariable <- V1 <- baseWeight <- calibWeight <- epsvalue <- f <- NULL
  temporary_hid <- temporary_hvar <- tmpVarForMultiplication <- value <- wValue <- wvst<- NULL
  dat_original <- dat
  dat <- copy(dat)
  nrowOriginal <- nrow(dat)
  ## originalsorting is fucked up without this
  dat[,OriginalSortingVariable:=.I]
  
  # dat sollte ein data.table sein
  # w ein Name eines Basisgewichts oder NULL
  ncp <- length(conP) # number of constraints on person level
  nch <- length(conH) # number of constraints on household level
  dimncp <- sapply(conP,function(x)prod(dim(x)))
  dimnch <- sapply(conH,function(x)prod(dim(x)))
  valueP <- paste0("valueP",seq_along(conP))###fixed target value, should not be changed in iterations
  valueH <- paste0("valueH",seq_along(conH))
  ###Housekeeping of the varNames used
  usedVarNames <- c(valueP,valueH,"value","baseWeight","wvst","wValue")
  renameVars <- NULL
  
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
  dat[, temporary_hid := as.factor(temporary_hid)]
  setnames(dat, "temporary_hid", hid)
  
  ## Names of the calibration variables for Person and household dimension
  pColNames <- lapply(conP,function(x)names(dimnames(x)))
  hColNames <- lapply(conH,function(x)names(dimnames(x)))
  
  for(i in seq_along(conP)){
    current_colnames <- pColNames[[i]]

    for(colname in current_colnames){
      if(!inherits(dat[[colname]], "factor")){
        message("converting column ", colname, " to factor.\n")
        set(
          dat, j = colname, 
          value = factor(dat[[colname]], levels = dimnames(conP[[i]])[[colname]])
        )
      }
      else if(!identical(levels(dat[[colname]]), dimnames(conP[[i]])[[colname]])){
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
        message("converting column ", colname, " to factor.\n")
        set(
          dat, j = colname, 
          value = factor(dat[[colname]], levels = dimnames(conH[[i]])[[colname]])
        )
      }
      else if(!identical(levels(dat[[colname]]), dimnames(conH[[i]])[[colname]])){
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
  
  pCalVar <- paste0("pcal",1:ncp)
  hCalVar <- paste0("hcal",1:nch)
  
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
        
        res <- calibP(i=i, conP=conP, epsP=epsP, dat=dat, error=error,
                      valueP=valueP, pColNames=pColNames,bound=bound, verbose=verbose, calIter=calIter, numericalWeighting=numericalWeighting)
        
        dat <- res[["dat"]]
        error <- res[["error"]]
        rm(res)
      }
      if(meanHH){
        ## replace person weight with household average
        dat[,calibWeight := geometric_mean(calibWeight, dat[[hid]])]
      }
      ### Household calib
      for(i in seq_along(conH)){
        res <- calibH(i=i, conH=conH, epsH=epsH, dat=dat, error=error,
                      valueH=valueH, hColNames=hColNames,bound=bound, verbose=verbose, calIter=calIter, looseH=looseH)
        dat <- res[["dat"]]
        error <- res[["error"]]
        rm(res)
      }
    }else{
      ### Person calib
      for(i in seq_along(conP)){
        
        res <- calibP(i=i, conP=conP, epsP=epsP, dat=dat, error=error,
                      valueP=valueP, pColNames=pColNames,bound=bound, verbose=verbose, calIter=calIter, numericalWeighting=numericalWeighting)
        dat <- res[["dat"]]
        error <- res[["error"]]
        rm(res)
        
        if(meanHH){
          ## replace person weight with household average
          dat[,calibWeight := geometric_mean(calibWeight, dat[[hid]])]
        }
        ### Household calib
        for(i in seq_along(conH)){
          res <- calibH(i=i, conH=conH, epsH=epsH, dat=dat, error=error,
                        valueH=valueH, hColNames=hColNames,bound=bound, verbose=verbose, calIter=calIter, looseH=looseH)
          dat <- res[["dat"]]
          error <- res[["error"]]
          rm(res)
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
  
  ## originalsorting is fucked up without this
  setkey(dat, OriginalSortingVariable)
  
  ## Return missings in calibWeight variable if no convergence was reached
  if(maxIter<calIter&returnNA){
    invisible(copy(dat_original)[, calibWeight := NA])
  }else{
    invisible(copy(dat_original)[,calibWeight := dat$calibWeight])  
  }  
}
#' @rdname ipu2
#' @export computeFrac
computeFrac <- function(curValue,Value,numericVar,weightVec){
  Value/curValue
}

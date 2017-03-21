#' Kish Factor
#' 
#' Compute the kish factor for a specific weight vector
#'
#' @name kishFactor
#' @param w a numeric vector with weights
#' @return The function will return the the kish factor 
#' @author Alexander Kowarik
#' @export
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
#' @return The function will return the input data \code{dat} with the
#' calibrated weights \code{calibWeight} as an additional column.
#' @seealso \code{\link{ipu}}
#' @export ipu2
#' @export computeLinear
#' @export computeFrac
#' @author Alexander Kowarik
#' @examples
#' data(eusilcS)
#' setDT(eusilcS)
#' eusilcS <- eusilcS[, list(db030,hsize,db040,age,rb090,netIncome,db090,rb050)]
#' 
#' ## some recoding
#' # generate age groups
#' eusilcS[age<0, age:=0]
#' eusilcS[,agegroup:=floor(age/10)]
#' # some recoding of netIncome for reasons of simplicity
#' eusilcS[is.na(netIncome), netIncome:=0] 
#' eusilcS[netIncome<0, netIncome:=0] 
#' # set hsize to 1,...,5+
#' eusilcS[hsize>=5, hsize:=5] 
#' 
#' ## example for base weights assuming a simple random sample of households stratified per region
#' eusilcS[, regSamp:=.N, by=db040]
#' eusilcS[, regPop:=sum(rb050), by=db040]
#' eusilcS[, baseWeight:=regPop/regSamp]
#' 
#' ## constraints on person level
#' # age 
#' conP1 <- xtabs(V1 ~ agegroup, data=eusilcS[,sum(rb050),by=agegroup])
#' # gender by region
#' conP2 <- xtabs(V1 ~ rb090+db040, data=eusilcS[,sum(rb050),by=list(rb090,db040)])
#' # personal net income by gender
#' conP3 <- xtabs(V1 ~ rb090, data=eusilcS[,sum(rb050*netIncome),by=rb090])
#' ## constraints on household level
#' conH1 <- xtabs(V1 ~ hsize+db040, data=eusilcS[!duplicated(db030),sum(rb050),list(hsize,db040)])
#' 
#' # array of convergence limits for conH1
#' epsH1 <- conH1
#' epsH1[as.character(1:4),] <- 0.005
#' epsH1["5",] <- 0.2
#' 
#' # without array epsP1
#' calibweights1 <- ipu2(eusilcS, hid = "db030", 
#'                       conP = list(conP1,conP2,netIncome=conP3), 
#'                       conH = list(conH1), 
#'                       epsP = list(1e-06,1e-06,1e-03),
#'                       epsH = 0.01,  
#'                       bound = NULL, verbose=TRUE,  maxIter = 200)
#' 
#' # with array epsP1, base weights and bound
#' calibweights2 <- ipu2(eusilcS, hid = "db030", 
#'                       conP = list(conP1,conP2), 
#'                       conH = list(conH1), 
#'                       epsP = 1e-06,
#'                       epsH = list(epsH1),  
#'                       w="baseWeight",
#'                       bound = 4, verbose=TRUE,  maxIter = 200)
#   fn <- function(a){
#    f <- a[1]*var+a[2]
#    (sum(f*var*w)-v)^2+(sum(f*w)-sum(w))^2
#  }
#  coef <- optim(c(1,1),fn)$par
ipu2 <- function(dat,hid=NULL,conP=NULL,conH=NULL,epsP=1e-6,epsH=1e-2,verbose=FALSE,
                 w=NULL,bound=4,maxIter=200,meanHH=TRUE,returnNA=TRUE,looseH=FALSE,numericalWeighting=computeLinear){
  
  OriginalSortingVariable <- V1 <- baseWeight <- calibWeight <- epsvalue <- f <- NULL
  temporary_hid <- temporary_hvar <- tmpVarForMultiplication <- value <- wValue <- wvst<- NULL
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
  
  mconP <- lapply(conP,melt,as.is=TRUE)##convert tables to long form
  mconH <- lapply(conH,melt,as.is=TRUE) 
  
  
  for(i in seq_along(conP)){
    # Harmonize the class of columns coming from the constraints
    # the result from melt is taken as character and then the class is set
    # melt does not make any sense for one dimensional numerical constraints however
    if(!(dim(conP[[i]])==1 && identical(colnames(mconP[[i]]),c("Var1", "value")))){
      cn <- colnames(mconP[[i]])[-ncol(mconP[[i]])]
      for(j in seq_along(cn)){
        cl <- class(dat[[cn[j]]])
        if("factor"%in%cl){
          mconP[[i]][[cn[j]]] <- factor(mconP[[i]][[cn[j]]],levels=levels(dat[[cn[j]]]))
        }else if("numeric"%in%cl){
          mconP[[i]][[cn[j]]] <- as.numeric(mconP[[i]][[cn[j]]])
        }else if("integer"%in%cl){
          mconP[[i]][[cn[j]]] <- as.integer(mconP[[i]][[cn[j]]])
        }
      }
      
      dat <- merge(dat,mconP[[i]],by=colnames(mconP[[i]])[-ncol(mconP[[i]])],all.x=TRUE,all.y=FALSE)
      setnames(dat,"value",valueP[i])
    }else{
      dat[,value:=conP[[i]]]
      setnames(dat,"value",valueP[i])
    }
  }
  for(i in seq_along(conH)){
    if(!(dim(conH[[i]])==1 && identical(colnames(mconH[[i]]),c("Var1", "value")))){
      cn <- colnames(mconH[[i]])[-ncol(mconH[[i]])]
      for(j in seq_along(cn)){
        cl <- class(dat[[cn[j]]])
        if("factor"%in%cl){
          mconH[[i]][[cn[j]]] <- factor(mconH[[i]][[cn[j]]],levels=levels(dat[[cn[j]]]))
        }else if("numeric"%in%cl){
          mconH[[i]][[cn[j]]] <- as.numeric(mconH[[i]][[cn[j]]])
        }else if("integer"%in%cl){
          mconH[[i]][[cn[j]]] <- as.integer(mconH[[i]][[cn[j]]])
        }
      }
      dat <- merge(dat,mconH[[i]],by=colnames(mconH[[i]])[-ncol(mconH[[i]])],all.x=TRUE,all.y=FALSE)
      setnames(dat,"value",valueH[i])
    }else{
      dat[,value:=conH[[i]]]
      setnames(dat,"value",valueH[i])
    }
  }
  
  if(nrow(dat)!=nrowOriginal){
    stop("There were problems merging the constraints to the data!\n")
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
  ## Names of the calibration variables for Person and household dimension
  pColNames <- lapply(conP,function(x)names(dimnames(x)))
  hColNames <- lapply(conH,function(x)names(dimnames(x)))
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
  
  if(is.list(epsP)){
    #mepsP <- lapply(epsP,melt,as.is=TRUE,value.name="epsvalue")##convert tables to long form
    mepsP <- lapply(epsP,meltepsfun)
    for(i in seq_along(epsP)){
      if(is.array(epsP[[i]])){
        cn <- colnames(mepsP[[i]])[-ncol(mepsP[[i]])]
        for(j in seq_along(cn)){
          cl <- class(dat[[cn[j]]])
          if("factor"%in%cl){
            mepsP[[i]][[cn[j]]] <- factor(mepsP[[i]][[cn[j]]],levels=levels(dat[[cn[j]]]))
          }else if("numeric"%in%cl){
            mepsP[[i]][[cn[j]]] <- as.numeric(mepsP[[i]][[cn[j]]])
          }else if("integer"%in%cl){
            mepsP[[i]][[cn[j]]] <- as.integer(mepsP[[i]][[cn[j]]])
          }
        }
      }
    }
  }
  if(is.list(epsH)){
    #mepsH <- lapply(epsH,melt,as.is=TRUE,value.name="epsvalue")##convert tables to long form
    mepsH <- lapply(epsH,meltepsfun)
    for(i in seq_along(epsH)){
      if(is.array(epsH[[i]])){
        cn <- colnames(mepsH[[i]])[-ncol(mepsH[[i]])]
        for(j in seq_along(cn)){
          cl <- class(dat[[cn[j]]])
          if("factor"%in%cl){
            mepsH[[i]][[cn[j]]] <- factor(mepsH[[i]][[cn[j]]],levels=levels(dat[[cn[j]]]))
          }else if("numeric"%in%cl){
            mepsH[[i]][[cn[j]]] <- as.numeric(mepsH[[i]][[cn[j]]])
          }else if("integer"%in%cl){
            mepsH[[i]][[cn[j]]] <- as.integer(mepsH[[i]][[cn[j]]])
          }
        }
      }
    }
  }
  
  ###Calib
  error <- TRUE
  calIter <- 1
  
  while(error&&calIter<=maxIter){
    error <- FALSE
    
    ### Person calib
    for(i in seq_along(conP)){
      if(is.list(epsP)){
        epsPcur <- epsP[[i]]
        mepsPcur <- mepsP[[i]]
      }else{
        epsPcur <- epsP
        # mepsPcur <- mepsP
      }
      
      if(isTRUE(names(conP)[i]!="")){
        ## numerical variable to be calibrated
        ## use name of conP list element to define numerical variable
        setnames(dat,names(conP)[i],"tmpVarForMultiplication")
        dat[,wValue:=sum(calibWeight*tmpVarForMultiplication),by=eval(pColNames[[i]])]
        setnames(dat,valueP[i],"value")
        # try to divide the weight between units with larger/smaller value in the numerical variable linear
        dat[,f:=numericalWeighting(head(wValue,1),head(value,1),tmpVarForMultiplication,calibWeight),by=eval(pColNames[[i]])]
        setnames(dat,"tmpVarForMultiplication",names(conP)[i])
        if(is.array(epsPcur)){## for numeric variables not the factor f is used but the abs relative deviation is computed per class
          curEps <- abs(dat[!is.na(f),1/(value/wValue)-1]) ## curEps is computed for all observations to compare it with the right epsValue
        }else{
          curEps <- dat[!is.na(f),max(abs(1/(value/wValue)-1))]  ## only the max curEps is needed because it is the same for all classed of the current variable
        }
      }else{
        # categorical variable to be calibrated
        dat[,wValue:=sum(calibWeight),by=eval(pColNames[[i]])]
        setnames(dat,valueP[i],"value")
        dat[,f:= value/wValue,by=eval(pColNames[[i]])]
        if(is.array(epsPcur)){
          curEps <- abs(dat[!is.na(f),1/f-1]) ## curEps is computed for all observations to compare it with the right epsValue
        }else{
          curEps <- dat[!is.na(f),max(abs(1/f-1))]  ## only the max curEps is needed because it is the same for all classed of the current variable
        }
      }
      if(is.array(epsPcur)){
        dat <- merge(dat,mepsPcur,by=pColNames[[i]],all.x=TRUE,all.y=FALSE)
        epsPcur <-dat[!is.na(f),epsvalue]
        dat[,epsvalue:=NULL]
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
      
      # if(any(curEps>epsPcur)){ ##  gleich zu oberer if-Abfrage dazugegeben
      #   error <- TRUE
      # }
      
      if(verbose&&any(curEps>epsPcur)&&calIter%%10==0){
        if(calIter%%100==0)
          print(subset(dat,!is.na(f))[abs(1/f-1)>epsPcur][,list(mean(f),.N),by=eval(pColNames[[i]])])
        #print(dat[abs(1/f-1)>epsPcur][,list(mean(f),.N),by=eval(pColNames[[i]])])
        cat(calIter, ":Not yet converged for P-Constraint",i,"\n")
      }
      
    }
    if(meanHH){
      dat[,calibWeight:=mean(calibWeight),by=eval(hid)] ## das machen wir bei MZ-HR-Paper vor der hh-Kalibrierung. Hier wird nur erstes hh-member kalibriert.
    }
    ### Household calib
    for(i in seq_along(conH)){
      if(is.list(epsH)){
        epsHcur <- epsH[[i]]
        mepsHcur <- mepsH[[i]]
      }else{
        epsHcur <- epsH
        #mepsHcur <- mepsH
      }
      
      if(isTRUE(names(conH)[i]!="")){
        ## numerical variable to be calibrated
        ## use name of conH list element to define numerical variable
        setnames(dat,names(conH)[i],"tmpVarForMultiplication")
        dat[,wValue:=sum(calibWeight*wvst*tmpVarForMultiplication),by=eval(hColNames[[i]])]
        setnames(dat,"tmpVarForMultiplication",names(conH)[i])
      }else{
        # categorical variable to be calibrated
        dat[,wValue:=sum(calibWeight*wvst),by=eval(hColNames[[i]])]
      }
      
      setnames(dat,valueH[i],"value")
      dat[,f:= value/wValue,by=eval(hColNames[[i]])]
      
      if(is.array(epsHcur)){
        dat <- merge(dat,mepsHcur,by=hColNames[[i]],all.x=TRUE,all.y=FALSE)
        curEps <- dat[,abs(1/f-1)]
        epsHcur <- dat[,epsvalue]
        dat[,epsvalue:=NULL]
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
            dat[,calibWeight:=boundsFakHH(g1=calibWeight,g0=baseWeight,eps=epsvalue,orig=value,p=wValue,bound=bound)]  
          }
        }else{
          dat[,calibWeight:=f*calibWeight,by=eval(hColNames[[i]])]
        }
      }
      # if("epsvalue"%in%colnames(dat)){
      #   dat[,epsvalue:=NULL]  
      # }
      
      setnames(dat,"value",valueH[i])
      if(any(curEps>epsHcur)){ 
        error <- TRUE
      }
      if(verbose&&any(curEps>epsHcur)&&calIter%%10==0){
        cat(calIter, ":Not yet converged for H-Constraint",i,"\n")
      }
    }
    if(verbose&&!error){
      cat("Convergence reached in ",calIter," steps \n")
    }else if(verbose&&maxIter==calIter){
      cat("Not converged in",maxIter,"steps \n")
    }
    calIter <- calIter + 1 
  }
  
  if(!is.null(w)){
    setnames(dat,"baseWeight",w)  
  }
  
  ##Housekeeping
  ###Housekeeping of the varNames used
  delVars <- c(delVars,"wvst","wValue","f")
  if(any(valueP%in%names(dat)))
    delVars <- c(delVars,valueP)
  if(any(valueH%in%names(dat)))
    delVars <- c(delVars,valueH)
  dat[,eval(delVars):=NULL]
  if(!is.null(renameVars)){
    setnames(dat,paste0(renameVars,"_safekeeping"),renameVars)
  }
  # loeschen; da macht man doch das umbenennen wieder rueckgaengig und das will man nicht!  
  #   if(any(names(dat)%in%usedVarNames)){
  #     renameVars <- names(dat)[names(dat)%in%usedVarNames]
  #     setnames(dat,renameVars,paste0(renameVars,"_safekeeping"))
  #   }
  ## originalsorting is fucked up without this
  setkey(dat,OriginalSortingVariable)
  dat[,OriginalSortingVariable:=NULL]
  
  ## Return missings in calibWeight variable if no convergence was reached
  if(maxIter<calIter&returnNA){
    invisible(dat[,calibWeight:=NA])  
  }else{
    invisible(dat)  
  }  
}
computeLinear <- function(curValue,Value,numericVar,weightVec,boundLinear=10){#current summed up value, correct summed up value, numeric variable, current weight
  h <- sum(weightVec*numericVar)
  j <- sum(weightVec*numericVar^2)
  N <- sum(weightVec)
  b <- (Value-N*j/h)/((-N*j/h)+h)
  a <- (N-b*N)/h
  f <- a*numericVar+b
  f[f<(1/boundLinear)] <- 1/boundLinear
  f[f>boundLinear] <- boundLinear
  return(f)
}
computeFrac <- function(curValue,Value,numericVar,weightVec){
  Value/curValue
}

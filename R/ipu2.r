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
#' Iterative Proportional Updating
#' 
#' Adjust sampling weights to given totals based on household-level and/or
#' individual level constraints.
#' 
#'
#' @name ipu2 
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
#' @return The function will return the input data \code{dat} with the
#' calibrated weights \code{calibWeight} as an additional column.
#' @seealso \code{\link{ipu}}
#' @export
#' @author Alexander Kowarik
#' @examples
#' nrowInp <- 1000
#' inp <- as.data.table(matrix(0, nrow=nrowInp, ncol=6))
#' colnames(inp) <- c("pid","hhid","income","agegroup","gender","region")
#' inp[,pid:=.I]
#' inp$hhid <- round(seq(from=1, to=round(nrowInp/3), length.out=nrowInp))
#' set.seed(223344)
#' # 4 age groups
#' inp$agegroup <- sample(c(1:4),nrowInp,replace=TRUE, prob=c(0.15,0.3,0.35,0.2))
#' # 2 gender groups
#' inp$gender <- sample(c(1,2),nrowInp,replace=TRUE)
#' # 3 regions
#' inp$region <- sample(c(1:3),nrowInp,replace=TRUE, prob=c(0.2,0.5,0.3))
#' inp[,region:=region[1], by=hhid]
#' # income
#' inp$income <- round(rpareto(nrowInp, 1000, 3))
#' inp[agegroup==1,income:=0]    
#' # hhincome
#' inp[,hhinc:=sum(income),by=list(hhid)]
#' 
#' # constraints on person level
#' conP1 <- array(c(239741,601386,360193,480699,1199962,718892,
#'   560069,1399490,840041,320257,799359,479911),c(3,4))
#' dimnames(conP1) <- list(region=c("1", "2", "3"),agegroup=c("1", "2", "3", "4"))
#' 
#' conP2 <- array(c(800171,1999596,1198754,800595,2000601,1200283),c(3,2))
#' dimnames(conP2) <- list(region=c("1", "2", "3"),gender=c("1", "2"))
#' 
#' conP3 <- array(c(1020162500,2549062073,1528896458,1021294346,2551562798,1530314970),c(3,2))
#' dimnames(conP3) <- list(region=c("1", "2", "3"),gender=c("1", "2"))
#' 
#' # constraints on household level
#' conH1 <- array(c(2041271296,5104297066,3055724783),3)
#' dimnames(conH1) <- list(region=c("1", "2", "3"))
#' 
#' conH2 <- array(c(533267,1333931,799469),3)
#' dimnames(conH2) <- list(region=c("1", "2", "3"))
#' 
#' # array of convergence limits for conP1 
#' epsP1 <- array(rep(c(0.9,rep(0.7,2)),4),c(3,4))
#' dimnames(epsP1) <- dimnames(conP1)
#' 
#' # example for base weights assuming a simple random sample of households stratified per region
#' inp[,xx:=1]
#' inp[, regSamp:=sum(xx),by=region]
#' reg <- cbind(region=1:3,regPop=addmargins(conP2,2)[1:3,"Sum"])
#' inp <- merge(inp,reg,by="region",sort=FALSE)
#' inp[,baseWeight:=regPop/regSamp]
#' inp[,xx:=NULL]
#' 
#' res1 <- ipu2(dat=inp, hid="hhid", w=NULL, conP=list(conP1, conP2, income=conP3),
#'   conH=list(hhinc=conH1, conH2), 
#'   epsP=0.09, epsH=0.05, verbose=TRUE, bound=NULL, maxIter=200, meanHH=TRUE)
#' 
#' # with array epsP1
#' res2 <- ipu2(dat=inp, hid="hhid", w=NULL, conP=list(conP1, conP2, income=conP3),
#'   conH=list(hhinc=conH1, conH2), 
#'   epsP=list(epsP1,0.07,0.07), epsH=0.05, verbose=TRUE, bound=NULL,
#'   maxIter=200, meanHH=TRUE)
#' 
#' # with base weights and bound
#' res3 <- ipu2(dat=inp, hid="hhid", w="baseWeight", conP=list(conP1, conP2),
#'   conH=list(hhinc=conH1, conH2), 
#'   epsP=list(epsP1,0.05), epsH=1e-2, verbose=TRUE, bound=4, maxIter=200,
#'   meanHH=TRUE)

ipu2 <- function(dat,hid=NULL,conP=NULL,conH=NULL,epsP=1e-6,epsH=1e-2,verbose=FALSE,
                 w=NULL,bound=4,maxIter=200,meanHH=TRUE,returnNA=TRUE){
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
  mconP <- lapply(conP,melt)##convert tables to long form
  mconH <- lapply(conH,melt)
  
  for(i in seq_along(conP)){
    dat <- merge(dat,mconP[[i]],by=colnames(mconP[[i]])[-ncol(mconP[[i]])])#,all.x=TRUE,all.y=FALSE)	
    setnames(dat,"value",valueP[i])
  }
  for(i in seq_along(conH)){
    dat <- merge(dat,mconH[[i]],by=colnames(mconH[[i]])[-ncol(mconH[[i]])])	
    setnames(dat,"value",valueH[i])
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
  ###Calib
  error <- TRUE
  calIter <- 1
  
  while(error&&calIter<=maxIter){
    error <- FALSE
    
    ### Person calib
    for(i in seq_along(conP)){
      if(is.list(epsP)){
        epsPcur <- epsP[[i]]
      }else{
        epsPcur <- epsP
      }
      
      if(isTRUE(names(conP)[i]!="")){
        ## numerical variable to be calibrated
        ## use name of conP list element to define numerical variable
        setnames(dat,names(conP)[i],"tmpVarForMultiplication")
        dat[,wValue:=sum(calibWeight*tmpVarForMultiplication),by=eval(pColNames[[i]])]    
        setnames(dat,"tmpVarForMultiplication",names(conP)[i])
      }else{
        # categorical variable to be calibrated
        dat[,wValue:=sum(calibWeight),by=eval(pColNames[[i]])] 
      }
      
      setnames(dat,valueP[i],"value")
      
      dat[,f:= value/wValue,by=eval(pColNames[[i]])] 
      
      if(is.array(epsPcur)){
        dat <- merge(dat,melt(epsPcur,value.name="epsvalue"),by=pColNames[[i]])
        curEps <- abs(dat[,1/f-1])
        epsPcur <- dat[,epsvalue]
        dat[,epsvalue:=NULL]
      }else{
        curEps <- dat[,max(abs(1/f-1))]
      }
      
      if(any(curEps>epsPcur)){## sicherheitshalber abs(epsPcur)? Aber es wird schon niemand negative eps Werte uebergeben??
        if(!is.null(bound)){
          dat[,calibWeight:=boundsFak(calibWeight,baseWeight,f,bound=bound)]#,by=eval(pColNames[[i]])]  
        }else{
          dat[,calibWeight:=f*calibWeight,by=eval(pColNames[[i]])]
        }
      }
      setnames(dat,"value",valueP[i])
      
      if(any(curEps>epsPcur)){ ## kann man doch eigentlich auch gleich zu oberer if-Abfrage dazugeben
        error <- TRUE
      }
      if(verbose&&(curEps>epsPcur)&&calIter%%10==0){
        if(calIter%%100==0)
          print(dat[abs(1/f-1)>epsPcur][,.(mean(f),.N),by=eval(pColNames[[i]])])
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
      }else{
        epsHcur <- epsH
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
        dat <- merge(dat,melt(epsHcur,value.name="epsvalue"),by=hColNames[[i]])
        curEps <- abs(dat[,1/f-1])
        epsHcur <- dat[,epsvalue]
        dat[,epsvalue:=NULL]
      }else{
        curEps <- dat[,max(abs(1/f-1))]
      }
      
      if(any(curEps>epsHcur)){    
        if(!is.null(bound)){
          dat[,calibWeight:=boundsFak(calibWeight,baseWeight,f,bound=bound),by=eval(hColNames[[i]])]  
        }else{
          dat[,calibWeight:=f*calibWeight,by=eval(hColNames[[i]])]
        }
      }
      setnames(dat,"value",valueH[i])
      
      if(any(curEps>epsHcur)){ 
      }
      if(verbose&&(curEps>epsHcur)&&calIter%%10==0){
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
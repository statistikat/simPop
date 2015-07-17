ipu <- function(inp, con, hid=NULL, eps=1e-07, verbose=FALSE) {
  if ( class(inp)[1] %in% c("data.frame", "data.table") ) {
    cnames <- colnames(inp)
    inp <- as.matrix(inp)
    colnames(inp) <- cnames
  }
  
  if ( is.null(hid) ) {
    hhid <- 1:nrow(inp)
  } else {
    ii <- match(hid, colnames(inp))
    hhid <- inp[,ii]
    inp <- inp[,-c(ii)]
  }
  inp <- inp[,match(names(con), colnames(inp))]
  con <- con[match(colnames(inp), names(con))]
  if ( any(colnames(inp) != names(con)) ) {
    stop("constraints (con) do not match input (inp) -> check your input!\n")
  }
  w <- as.numeric(rep(1,nrow(inp)))
  con <- as.numeric(unlist(con))
  w <- .Call("simPop_ipu_work", inp, con, w, eps=eps, verbose=ifelse(verbose, 1L, 0L), package="simPop")
  out <- cbind(hhid, inp)
  out <- as.data.frame(out)
  out$weights <- w
  invisible(out)
}
ipu2 <- function(dat,conP,conH,hid=NULL,epsP=1e-2,epsH=1e-2,verbose=FALSE,
    w=NULL,bound=4,maxIter=200,meanHH=TRUE){
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
  }
  ### Treatment of HID, creating 0,1 var for being the first hh member
  delVars <- c()
  if(is.null(hid)){
    delVars <- c("hid")
    hid <- "hid"
    dat[,hid:=1:nrow(dat)]
    dat[,wvst:=1] 
  }else{
    dat[,wvst:=as.numeric(!duplicated(dat[,hid,with=FALSE]))]  
  }
  
  
  boundsFak <- function(g1,g0,f,bound=4){ # Berechnet die neuen Gewichte (innerhalb 4, .25 Veraenderungsraten)
    g1 <- g1 * f
    g1[(g1/g0)>bound] <- bound*g0[(g1/g0)>bound]
    g1[(g1/g0)<(1/bound)] <- (1/bound)*g0[(g1/g0)<(1/bound)]
    return(g1)
  }
  mconP <- lapply(conP,melt)##convert tables to long form
  mconH <- lapply(conH,melt)
  
  for(i in seq_along(conP)){
    dat <- merge(dat,mconP[[i]],by=colnames(mconP[[i]])[-ncol(mconP[[i]])])	
    setnames(dat,"value",valueP[i])
  }
  for(i in seq_along(conH)){
    dat <- merge(dat,mconH[[i]],by=colnames(mconH[[i]])[-ncol(mconH[[i]])])	
    setnames(dat,"value",valueH[i])
  }
  pCalVar <- paste0("pcal",1:ncp)
  hCalVar <- paste0("hcal",1:nch)
  if(is.null(w)){
    if(!is.null(bound)&&is.null(w))
      stop("Bounds are only reasonable if base weights are provided")
    dat[,calibWeight:=1]
    delVars <- c(delVars,"baseWeight")
  }else{
    dat[,calibWeight:=dat[,w,with=FALSE]]
    setnames(dat,w,"baseWeight")
  }
  ## Names of the calibration variables for Person and household dimension
  pColNames <- lapply(conP,function(x)names(dimnames(x)))
  hColNames <- lapply(conH,function(x)names(dimnames(x)))
  ###Calib
  error <- TRUE
  calIter <- 1
  while(error&&calIter<=maxIter){
    error <- FALSE
    
    ### Person calib
    for(i in seq_along(conP)){
      dat[,wValue:=sum(calibWeight),by=eval(pColNames[[i]])]
      setnames(dat,valueP[i],"value")
      dat[,f:=value/wValue,by=eval(pColNames[[i]])]
      if(!is.null(bound)){
        dat[,calibWeight:=boundsFak(calibWeight,baseWeight,f,bound=bound),by=eval(pColNames[[i]])]  
      }else{
        dat[,calibWeight:=f*calibWeight,by=eval(pColNames[[i]])]
      }
      setnames(dat,"value",valueP[i])
      curEps <- dat[,max(abs(f-1))]
      if(curEps>epsP){
        error <- TRUE
      }
      if(verbose&&curEps&&calIter%%10==0){
        cat(calIter, ":Not yet converged for P-Constraint",i,"\n")
      }
    }
    ### Household calib
    for(i in seq_along(conH)){
      dat[,wValue:=sum(calibWeight*wvst),by=eval(hColNames[[i]])]
      setnames(dat,valueH[i],"value")
      dat[,f:=value/wValue,by=eval(hColNames[[i]])]
      if(!is.null(bound)){
        dat[,calibWeight:=boundsFak(calibWeight,baseWeight,f,bound=bound),by=eval(hColNames[[i]])]  
      }else{
        dat[,calibWeight:=f*calibWeight,by=eval(hColNames[[i]])]
      }
      setnames(dat,"value",valueH[i])
      curEps <- dat[,max(abs(f-1))]
      if(curEps>epsH){
        error <- TRUE
      }
      if(verbose&&curEps&&calIter%%10==0){
        cat(calIter, ":Not yet converged for H-Constraint",i,"\n")
      }
      if(meanHH){
        dat[,calibWeight:=mean(calibWeight),by=eval(hid)]
      }
    }
    calIter <- calIter + 1
    if(verbose&&!error){
      cat("Convergence reached in ",calIter," steps \n")
    }else if(verbose&&maxIter==calIter){
      cat("Not converged in",maxIter,"steps \n")
    }
  }
  setnames(dat,"baseWeight",w)
  ##Housekeeping
  ###Housekeeping of the varNames used
  delVars <- c(delVars,valueP,valueH,"wvst","wValue","f")
  dat[,eval(delVars):=NULL]
  if(!is.null(renameVars)){
    setnames(dat,paste0(renameVars,"_safekeeping"),renameVars)
  }
  if(any(names(dat)%in%usedVarNames)){
    renameVars <- names(dat)[names(dat)%in%usedVarNames]
    setnames(dat,renameVars,paste0(renameVars,"_safekeeping"))
  }
  invisible(dat)  
}



#rm(list=ls())
#source("/Users/alex/git/simpopulation2/R/ipu.R")
#require(simPop);require(reshape2)
#data(eusilcS)
#eusilcS$agecut <- cut(eusilcS$age, 7)
#eusilcS$nat <- sample(LETTERS[1:5],nrow(eusilcS),replace=TRUE)
#eusilcS$emp <- sample(LETTERS[1:5],nrow(eusilcS),replace=TRUE)
#totals1 <- tableWt(eusilcS[, c("rb090","agecut","db040")], weights=eusilcS$rb050)
#totals2 <- tableWt(eusilcS[, c("rb090","emp","db040")], weights=eusilcS$rb050)
#totals3 <- tableWt(eusilcS[, c("nat","db040")], weights=eusilcS$rb050)
#totals1h <- tableWt(eusilcS[!duplicated(eusilcS$db030), c("hsize","db040")], weights=eusilcS$rb050[!duplicated(eusilcS$db030)])
#conP <- list(totals1,totals2,totals3)
#conH <- list(totals1h)
#dat <- data.table(eusilcS)
#cal1 <- ipu2(dat,conP=list(totals1,totals2,totals3),conH=list(totals1h),verbose=TRUE,
#    w="rb050",bound=4,maxIter=200,meanHH=TRUE,hid="db030")
#
#cal2 <- ipu2(dat,conP=list(totals1,totals2,totals3),conH=list(totals1h),verbose=TRUE,
#    w="rb050",bound=4,maxIter=200,meanHH=FALSE,hid="db030")

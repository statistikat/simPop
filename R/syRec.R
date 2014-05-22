syRec <- function(x, basisvar, hid, newvar, weights,
                  size=sum(weights)/nrow(x), strata){
  if(!is.factor(x[,newvar])) {
    stop("newvar is not a factor. synthetic reconstruction works only with factors.")
  }
  pop <- sampHH(x, sizefactor=size, 
                      hid="db030", strata="db040", 
                      hsize="hsize")[,c(basisvar,hid),with=FALSE]
  ## simulate new variable  
  stat <- as.numeric(levels(eusilcS[,newvar]))
  ## get conditional probablities:
  vars <- c(basisvar,newvar)
  TAB1 <- prop.table(table(eusilcS[,vars]))
  ## simulate new variable:
#  pop <- data.frame(age=eusilcS[eusilcS$db030==hh1,"age"],
#                    sex=eusilcS[eusilcS$db030==hh1,"rb090"])
  samp <- function(x){
    if(all(TAB1[as.numeric(x[1]),as.numeric(x[2]),]==0)){
      status <- NA
    } else {
      status <- sample(stat, 1, prob=TAB1[as.numeric(x[1]),as.numeric(x[2]),])
    }
    as.numeric(status)
  }
  pop <- data.frame(pop)
  stat <- as.numeric(levels(eusilcS[,newvar]))
  pop$status <- aaply(pop, 1, .expand=FALSE, .fun=function(x){
    if(all(TAB1[as.numeric(x[1]),as.numeric(x[2]),]==0)){
      status <- NA
    } else {
      status <- sample(stat, 1, 
                       prob=TAB1[as.numeric(x[1]),as.numeric(x[2]),])
    }
    as.numeric(status)
  })
  pop   
}


syRec <- function(synthPopObj, additional,
                  method=c("direct", "multinom", "distribution"), seed=1){
#  pop <- sampHH(x, sizefactor=size, 
#                hid="db030", strata="db040", 
#                hsize="hsize")[,c(basisvar,hid),with=FALSE]
  x <- NULL
  
  dataP <- synthPopObj@pop
  dataS <- synthPopObj@sample
  data_pop <- dataP@data
  data_sample <- dataS@data
  basic <- synthPopObj@basicHHvars
  
  if ( any(additional %in% colnames(data_pop)) ) {
    stop("variables already exist in the population!\n")
  }
  
  parallel <- FALSE
  have_win <- Sys.info()["sysname"] == "Windows"
  nr_cores <- detectCores()
  if ( nr_cores > 2 ) {
    parallel <- TRUE
    nr_cores <- nr_cores-1 # keep one core available
  } else {
    parallel <- FALSE
  }
  
  ##### initializations
  if ( !missing(seed) ) {
    set.seed(seed)  # set seed of random number generator
  }
  method <- match.arg(method)
  
  if ( method == "multinom" ) {
    basic <- c(basic, dataP@hhsize)
  }
  if ( !all(basic %in% colnames(data_pop)) ) {
    stop("undefined variables in the population data -> check your input!\n")
  }
  
  # check sample data against additional parameter
  if ( !all(additional %in% colnames(data_sample)) )  {
    stop("undefined variables in the sample data -> check your input!\n")
  }
  # variables are coerced to factors
  data_sample <- checkFactor(data_sample, c(dataS@strata, basic, additional))
  data_pop <- checkFactor(data_pop, c(dataP@strata, basic))
  
  # check arguments to account for structural zeros
  if ( length(additional) == 1 ) {
    if ( !(length(limit) == 1 && isTRUE(names(limit) == additional)) ) {
      limit <- list(limit)
      names(limit) <- additional
    }
    if ( !(length(censor) == 1 && isTRUE(names(censor) == additional)) ) {
      censor <- list(censor)
      names(censor) <- additional
    }
  }
  
  
  ## simulate new variable  
  stat <- levels(dataS[,additional]))
  ## get conditional probablities:
  vars <- c(basisvar,newvar)
  TAB1 <- prop.table(table(eusilcS[,vars]))
  ## simulate new variable:
  #  pop <- data.frame(age=eusilcS[eusilcS$db030==hh1,"age"],
  #                    sex=eusilcS[eusilcS$db030==hh1,"rb090"])
  samp <- function(x){
    if(all(TAB1[as.numeric(x[1]),as.numeric(x[2]),]==0)){
      status <- NA
    } else {
      status <- sample(stat, 1, prob=TAB1[as.numeric(x[1]),as.numeric(x[2]),])
    }
    as.numeric(status)
  }
  pop <- data.frame(pop)
  stat <- as.numeric(levels(eusilcS[,newvar]))
  pop$status <- aaply(pop, 1, .expand=FALSE, .fun=function(x){
    if(all(TAB1[as.numeric(x[1]),as.numeric(x[2]),]==0)){
      status <- NA
    } else {
      status <- sample(stat, 1, 
                       prob=TAB1[as.numeric(x[1]),as.numeric(x[2]),])
    }
    as.numeric(status)
  })
  pop   
}




install_bitbucket("simPopulation2", username="bernhard_da",  auth_user
                  = "matthias-da", password = "Hdjexly1")

library(synthPop)
require(plyr)
data(eusilcS)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize",
                    strata="db040", weight="db090")
synthPopObj <- simStructure(data=inp, method="direct",
                            basicHHvars=c("age", "rb090", "hsize", "pl030", "pb220a"))

weights <- eusilcS$rb050
selectHHs <- sampHH(eusilcS, sizefactor=sum(weights)/nrow(eusilcS), 
                         hid="db030", strata="db040", hsize="hsize")

syRec(eusilcS, basisvar=c("age","rb090"), hid="db030",
      newvar="pl030", weights="rb050", size=0.01, strata="db040")



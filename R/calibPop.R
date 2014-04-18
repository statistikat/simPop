obj <- function(data,totals,weights,parameter) {  
  index <- weights == 1
  synthtotals <- data[index,parameter, with=F][,.N, by=key(data)]
  setkeyv(synthtotals, parameter)
  setkeyv(totals, parameter)
  dat <- totals[synthtotals]
  return(sum(abs(dat$N - dat$N.1)))
}

## function for sample weights
## use difference to marginals
## sample weight equals positive difference between marginals and 
## synth-data divided by whole deviation
sample.weights <- function(data,totals,weights,parameter, donor=TRUE) {  
  table.weights <- sweight.donor <- sweight.synth <- NULL
  index <- weights == 1
  synthtotals <- data[index,parameter, with=F][,.N, by=key(data)]
  setkeyv(synthtotals, parameter)
  setkeyv(totals, parameter)
  dat <- totals[synthtotals]
  if ( !donor ) {
    dat[,diff:=dat$N.1-dat$N]
  } else {
    dat[,diff:=dat$N-dat$N.1]
  }
  
  ## let categories with too many ind be selected by a vary small chance
  ## use abs-total-diff^-1 as testing option
  pindex <- dat$diff<0
  dat$diff[pindex] <- sum(abs(dat$diff))^-1  
  dat[,table.weights:= dat$diff/sum(dat$diff)]
  dat[,c(split,"N","N.1","diff"):=NULL]
  
  #gibt ein WARINING falls table.weights noch nicht als col definiert
  if ( donor == TRUE ) {
    if ( !is.na(match("sweight.donor", colnames(data))) ) {
      data[,sweight.donor:=NULL]  
    }    
    data <- merge(data,dat,suffixes=c())
    setnames(data,"table.weights","sweight.donor")
  } else {
    if ( !is.na(match("sweight.synth", colnames(data))) ) {
      data[,sweight.synth:=NULL]
    }    
    data <- merge(data,dat,suffixes=c())
    setnames(data,"table.weights","sweight.synth")
  }  
  return(data)
}

## resample function
## switches a percentage of weights, given by factor, from 0 to 1
## switches full households always
resample <- function(data, factor) {
  data$id <- 1:nrow(data)
  new.weights <- data$weights
  active <- new.weights == 1  
  inactive <- new.weights == 0 
  
  randomactive <- sample(data$hid[active],factor,prob=data$sweight.synth[active])
  a <- data.table(hid=randomactive)
  setkeyv(a, "hid")
  activeindex <- data[a]$id
  new.weights[activeindex] <- 0   
  
  randominactive <- sample(data$hid[inactive],factor,prob=data$sweight.donor[inactive])  
  b <- data.table(hid=randominactive)
  setkeyv(b, "hid")  
  inactiveindex <- data[b]$id
  new.weights[inactiveindex] <- 1
  return(new.weights)
}

## data...data containing the households ect
## totals...marginals to be reached using simulated annealing
## temp...startin temperatur (use 30 for testing purpose)
## eps...acceptance level for obejctive function (use 500 for testing purpose)
## maxiter...max number of iterations during each temperature setting
## split...problem will be split into sub-problems based of the factorization of split
## parameter...vector which holds the attributs to be considered in the sa-algorithm   

# main workhorse, function because easier for parallelization
one_run <- function(data0, totals0, parameter, temp, eps.factor, 
  maxiter, temp.cooldown, factor.cooldown, min.temp, sample) {
  
  temp0 <- temp
  ##count number of temp coolings - for testing purpose
  cooldown <- 0 
  ## c count if sampleweights are swapped
  change <- 0
  
  setkey(data0, "hid")
  
  ## initialize weights-vector in accordance to small problem
  index0 <- data0[round(sum(totals0$N)),]$hid
  index1 <- max(which(c(data0$hid==index0)==TRUE))
  
  weights <- c(rep(1, index1),rep(0,nrow(data0)-index1))
  data0$weights <- weights
  data0$new.weights <- weights
  
  ## acceptable error
  eps <- eps.factor*sum(data0$weights)
  
  ## initial number of swaps per resample
  ## decreases with temperature
  factor <- index1*(1/10)  
  
  setkeyv(data0, parameter)
  objective <- obj(data=data0, totals=totals0, weights=data0$weights, parameter=parameter)
  
  if ( sample ) {
    data0 <- sample.weights(data=data0, totals=totals0, weights=data0$weights, parameter=parameter, donor=FALSE)
    data0 <- sample.weights(data=data0, totals=totals0, weights=data0$weights, parameter=parameter, donor=TRUE)
  } else {
    data0$sweight.synth <- c(rep(1,nrow(data0)))
    data0$sweight.donor <- c(rep(1,nrow(data0)))
  }
  
  if ( objective <= eps ) { 
    out <- list()
    text <- paste("synthetic data fullfils marginals to the error of", eps)    
    out$solution.list <- text
    out$objective.list <- c(objective,eps)
    out$weights.list <- data0$weights
    out$s.weights.list <- cbind(data0$sweight.synth,data0$sweight.donor)
    out$factor.list <- factor
    out$temp.list <- temp0
    out$cooldown.list <- cooldown
    out$changes.list <- change   
  } else {
    ## minimal temperature 
    ## if temperature falls below minimal temp -> terminate
    while( temp0 > min.temp ) {      
      n <- 1
      while( n<maxiter ) {    
        ## k count if resample does not give a better solution
        ## swap zeros to one in weights.new
        setkey(data0,"hid")
        data0$new.weights <- resample(data=data0, factor=factor)
        setkeyv(data0,parameter)
        objective.new <- obj(data=data0, totals=totals0, weights=data0$new.weights, parameter = parameter)
        
        ## if new sample fullfils marginals -> terminate
        if ( objective.new <= eps ) {
          objective <- objective.new
          data0$weights <- data0$new.weights
          change <- change+1
          break
        }
        
        ## choose wether to accepts the resample
        if ( objective.new <= objective ) { 
          data0$weights <- data0$new.weights 
          objective <- objective.new
          if ( sample ) {
            data0 <- sample.weights(data = data0, totals = totals0 , weights = data0$weights , parameter = parameter, donor=FALSE)
            data0 <- sample.weights(data = data0, totals = totals0 , weights = data0$weights , parameter = parameter, donor=TRUE)
          }
          change <- change+1
        }
        
        if ( objective.new > objective ) {      
          prob <- exp(-(objective.new-objective)/temp)
          x <- sample(c(0,1), 1,prob=c(1-prob,prob))
          
          if ( x == 1 ) { 
            data0$weights <- data0$new.weights 
            objective <- objective.new
            if ( sample ) {
              data0 <- sample.weights(data=data0, totals=totals0, weights=data0$weights, parameter=parameter, donor=FALSE)
              data0 <- sample.weights(data=data0, totals=totals0, weights=data0$weights, parameter=parameter, donor=TRUE)
            }
            change <- change+1
          }
        }    
        n <- n+1        
      }
      ## decrease temp and decrease factor accordingly
      ## decrease temp by a const fraction (simple method used for testing only)
      temp0 <- temp.cooldown*temp0
      factor <- floor(factor.cooldown*factor)
      if ( factor == 0 ) {
        factor <- 1
      }
      cooldown <- cooldown + 1
      if ( objective.new <= eps | cooldown == 500 ) {
        break
      }       
    }
    
    if ( objective <= eps ) {
      text <- paste("synthetic data fullfils marginals to the error of", eps,"- synthetic data indicated by weights")
    } else {
      text <- paste("simulated annealing failed to optimize objective function to the error of", eps,"- synthetic data indicated by weights")  
    }
    
    out <- list()
    out$solution.list <- text
    out$objective.list <- c(objective,eps)
    out$weights.list <- data0$weights
    out$s.weights.list <- cbind(data0$sweight.synth,data0$sweight.donor)
    out$factor.list <- factor
    out$temp.list <- temp0
    out$cooldown.list <- cooldown
    out$changes.list <- change    
  }  
  return(out)
}


###########################################################
###########################################################
## Simulated Annealing Function
## 
## A Simulated Annealing Algorithm for the Creation of Synthetic Microdata. Aims to find, given a microdataset, a combination of different households which optimaly satisfy, in the sense of an acceptable error, a given table of specific attributes.
##
##
## data = whole dataset
## totals = target-marginals
## hid = household-ID
## parameter = vector containg the attributes to be optimized for.
## split = vector containing attribute in which the Problem will be splitted.
## temp = starting temperature for the simulated annealing algorithm.
## eps.factor = scalar; percentage of populatin of synthetic-data which corresponds abort criterion.
## maxiter = number of iterations within a temperature step.
## temp.cooldown = scalar between 0 and 1; determines the cooldown step for the temperatur.
## factor.cooldown = scalar between 0 and 1; determines the cooldown step for the factor(resample.R).
## min.temp = minimal acceptable temperature; if temp < min.temp stop algorithm.
## sample = boolean variable; If TRUE sample weights will be used for resample algorithm. 
calibPop <- function(data, totals, hid, parameter, 
  split, temp = 30, eps.factor = 0.05, maxiter=200, temp.cooldown = 0.975, 
  factor.cooldown = 0.85, min.temp = 10^-3, sample = FALSE, parallel=FALSE) {
  
  table.weights <- sweight.donor <- sweight.synth <- NULL
  
  if ( parallel ) {
    if ( Sys.info()["sysname"] != "Windows" ) {
      nr_cores <- detectCores()
      if ( nr_cores > 2 ) {
        parallel <- TRUE
        nr_cores <- nr_cores-1 # keep one core available
      }    
    } else {
      parallel <- FALSE 
    }    
  }
  
  # calculate hhsize if not existing
  setkeyv(data, hid)
  res <- data[,.N,by=key(data)]
  setkeyv(res, hid)
  setnames(res, c(hid, "hhsize_calculated"))
  data <- res[data]
  
  ii <- match(parameter, colnames(data))
  data[,ii] <- data[,lapply(.SD,as.factor),.SDcols=parameter]
  rm(ii)
  
  ## split problem by "split"-factor
  split.number <- unique(data[,split,with=F])
  setkeyv(data,split) 
  
  ## initialize for output for every smaller problem    
  solution.list <- c()
  objective.list <- c()
  weights.list <- list()
  s.weights.list <- list()
  factor.list <- c()
  temp.list <- c()
  cooldown.list <- c()
  changes.list <- c()
  
  ## use sa-algorithm for every smaller problem individually  
  if ( parallel ) {
    results <- mclapply(1:nrow(split.number), function(x) { 
      one_run(
        data0=data[split.number[x]], 
        totals0=totals[which(totals[,split,with=F]==as.character(split.number[x][[split]])),], 
        parameter, temp, eps.factor, maxiter, 
        temp.cooldown, factor.cooldown, 
        min.temp, sample)
      }, 
      mc.cores = max(nr_cores,nrow(split.number)))    
    
    for ( i in 1:nrow(split.number) ) {
      solution.list <- c(solution.list, results[[i]]$solution.list)
      objective.list <- rbind(objective.list, results[[i]]$objective.list)  
      weights.list[[i]] <- results[[i]]$weights.list
      s.weights.list[[i]] <- results[[i]]$s.weights.list      
      factor.list <- c(factor.list, results[[i]]$factor.list)
      temp.list <- c(temp.list, results[[i]]$temp.list)
      cooldown.list <- c(cooldown.list, results[[i]]$cooldown.list)
      changes.list <- c(changes.list, results[[i]]$changes.list)          
    }    
  } else {
    for( i in 1:nrow(split.number) ) {    
      cat("i=",i,"\n")
      data0 <- data[split.number[i]]
      totals0 <- totals[which(totals[,split,with=F]==as.character(split.number[i][[split]])),]
      
      res <- one_run(data0, totals0, parameter, temp, eps.factor, maxiter, temp.cooldown, factor.cooldown, min.temp, sample)
      
      solution.list <- c(solution.list, res$solution.list)
      objective.list <- rbind(objective.list, res$objective.list)
      weights.list[[i]] <- res$weights.list
      s.weights.list[[i]] <- res$s.weights.list
      factor.list <- c(factor.list, res$factor.list)
      temp.list <- c(temp.list, res$temp.list)
      cooldown.list <- c(cooldown.list, res$cooldown.list)
      changes.list <- c(changes.list, res$changes.list)    
    }      
  }
  solution.list <- cbind(solution.list, split.number)
  objective.list <- cbind(objective.list, split.number)
  factor.list <- cbind(factor.list, split.number)
  temp.list <- cbind(temp.list, split.number)
  cooldown.list <- cbind(cooldown.list, split.number)
  changes.list <- cbind(changes.list, split.number)     
  result <- list(Solution = solution.list , Obj.value = objective.list, sample.weight = s.weights.list, weights = weights.list, parameter = parameter, factor = factor.list, factor.cooldown = factor.cooldown, cooldown = cooldown.list, temp.cooldown = temp.cooldown, Temperature = temp.list , number.changes = changes.list)
  return(invisible(result))
}

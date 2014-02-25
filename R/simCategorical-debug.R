# cur.var: current.variable (value of 'i' from outer loop)
generate.values <- function(dataSample, dataPop, method, cur.var, excludeLevels, hasNewLevels, newLevels, w, formula.cmd, eps, limit, censor, levelsResponse) {
  if( !nrow(dataSample) ) {
    return(character())
  }

  indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
  grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]

  # in sample, observations with NAs have been removed to fit the 
  # model, hence population can have additional levels
  # these need to be removed since those probabilities cannot 
  # be predicted from the model
  if( excludeLevels ) {
    exclude <- mapply(function(pop, new) pop %in% new, 
      pop=grid[, hasNewLevels, drop=FALSE], 
      new=newLevels[hasNewLevels])
    if(is.null(dim(exclude))) {
      exclude <- which(any(exclude))
    } else exclude <- which(apply(exclude, 1, any))
  } else exclude <- integer()

  # fit multinomial model
  # command needs to be constructed as string
  # this is actually a pretty ugly way of fitting the model
  #browser()
  tt <- Sys.time()
  mod <- eval(parse(text=formula.cmd))  # fitted model
  tt <- Sys.time()-tt
  print(tt)

  # rpart
  #tt <- Sys.time()
  #mod <- rpart(pl030 ~ age + rb090 + hsize, weights=rb050, data=dataSample)
  #tt <- Sys.time()-tt
  #print(tt)	

  # predict probabilities
  if(length(exclude) == 0) {
    if ( method %in% c("multinom", "rpart", "naivebayes") ) {
      probs <- predict(mod, newdata=grid, type="probs")
    }
    if ( method == "liblinear" ) {
      probs <- predict(mod, newx=grid, proba=TRUE)
    }
  } else {
    if ( method %in% c("multinom", "rpart", "naivebayes") ) {
      probs <- predict(mod, newdata=grid[-exclude, , drop=FALSE], type="probs")
    }
    if ( method == "liblinear" ) {
      probs <- predict(mod, newx=grid[-exclude, , drop=FALSE], proba=TRUE)
    }
  }
  # set too small probabilities to exactly 0
  if( !is.null(eps) ) {
    probs[probs < eps] <- 0
  }

  # ensure code works for missing levels of response
  ind <- as.integer(which(table(dataSample[, cur.var]) > 0))
  if( length(ind) > 2 && (nrow(grid)-length(exclude)) == 1 ) {
    probs <- t(probs)
  }

  # account for structural zeros
  if ( (!is.null(limit) || !is.null(censor)) && !is.null(dim(probs)) ) {
    if(length(exclude) == 0) {
      probs <- adjustProbs(probs, grid, names(indGrid), limit[[i]], censor[[i]])
    } else {
      probs <- adjustProbs(probs, grid[-exclude, , drop=FALSE], names(indGrid)[-exclude], limit[[i]], censor[[i]])
    }
  }
  # local function for sampling from probabilities
  if( length(ind) == 1 ) {
    resample <- function(k, n, p) rep.int(1, n[k])
  } else if(length(ind) == 2) {
    resample <- function(k, n, p) spSample(n[k], c(1-p[k],p[k]))
  } else resample <- function(k, n, p) spSample(n[k], p[k,])
  # generate realizations for each combination
  if(length(exclude) == 0) {
    ncomb <- as.integer(sapply(indGrid, length))
    sim <- lapply(1:length(ncomb), resample, ncomb, probs)
  } else {
    ncomb <- as.integer(sapply(indGrid[-exclude], length))
    sim <- as.list(rep.int(NA, length(indGrid)))
    sim[-exclude] <- lapply(1:length(ncomb), resample, ncomb, probs)
  }
  sim <- unsplit(sim, dataPop, drop=TRUE)
  # return realizations
  levelsResponse[ind][sim]
}

#x <- "Burgenland"
#i <- "pl030"

#l <- unique(dataS$db040)
#tt <- Sys.time()
#for ( x in l ) {
#	system.time(out <- generate.values(
#		dataSample=dataS[dataS[, strata] == x, , drop=FALSE], 
#		dataPop=dataP[indStrata[[x]], predNames, drop=FALSE],
#		cur.var=i,
#		excludeLevels=excludeLevels, 
#		hasNewLevels=hasNewLevels, 
#		newLevels=newLevels,
#		w=w,formula.cmd=formula.cmd,
#		eps=eps,limit=limit,censor=censor,
#		levelsResponse=levelsResponse))	
#}
#Sys.time() - tt

simCategorical2 <- function(dataS, dataP, method, parallel=TRUE) {
  ##### initializations
  w = "rb050"
  strata = "db040"
  basic <- c("age", "rb090", "hsize")
  additional = c("pl030", "pb220a")
  method = "multinom"
  limit <- NULL
  censor <- NULL
  maxit = 500
  MaxNWts = 1500
  eps = NULL
  seed <- 10

  varNames <- c(w=w, strata=strata, basic, additional)

  time.prep <- Sys.time()

  # check data
  if(all(varNames %in% names(dataS))) {
    dataS <- dataS[, varNames]
  } else stop("undefined variables in the sample data")
  if(!all(c(strata, basic) %in% names(dataP))) {
    stop("undefined variables in the population data")
  }

  # observations with missings are excluded from simulation
  exclude <- getExclude(dataS)
  if(length(exclude)) dataS <- dataS[-exclude,]

  # variables are coerced to factors
  dataS <- checkFactor(dataS, c(strata, basic, additional))
  dataP <- checkFactor(dataP, c(strata, basic))

  # check arguments to account for structural zeros
  if(length(additional) == 1) {
    if(!(length(limit) == 1 && isTRUE(names(limit) == additional))) {
      limit <- list(limit)
      names(limit) <- additional
    }
    if(!(length(censor) == 1 && isTRUE(names(censor) == additional))) {
      censor <- list(censor)
      names(censor) <- additional
    }
  }

  # list indStrata contains the indices of dataP split by strata
  N <- nrow(dataP)
  indStrata <- split(1:N, dataP[, strata, drop=FALSE])

  time.prep <- Sys.time() - time.prep
  cat("time.prep: "); print(time.prep)

  ##### simulation
  # predictor variables
  predNames <- basic  # names of predictor variables

  #cat("additional:", additional, "\n")
  for(i in additional ) {
    cat("i:",i,"\n")
    # components of multinomial model are specified
    levelsResponse <- levels(dataS[, i])

    # simulation of variables using a sequence of multinomial models
    if ( method == "multinom" ) {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("suppressWarnings(multinom(", formula.cmd, 
        ", weights=", w, ", data=dataSample, trace=FALSE", 
        ", maxit=",maxit, ", MaxNWts=", MaxNWts,"))", sep="")					
    }
    # simulation via recursive partitioning and regression trees
    if ( method == "rpart" ) {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("suppressWarnings(rpart(", formula.cmd, ", weights=", w, ", data=dataSample, method='class'))", sep="")						
    }
    if ( method == "naivebayes" ) {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("suppressWarnings(NaiveBayes(", formula.cmd, ", data=dataSample, usekernel=TRUE))", sep="")						
    }	
    if ( method == "liblinear" ) {
      formula.cmd <- paste("suppressWarnings(LiblineaR(data=dataS[,predNames], labels=dataS[,i], type=0, bias=TRUE,verbose=FALSE))", collapse="")						
    }				

    time.levs <- Sys.time()

    # check if population data contains factor levels that do not exist 
    # in the sample
    newLevels <- lapply(predNames, function(nam) {
      levelsS <- levels(dataS[, nam])
      levelsP <- levels(dataP[, nam])
      levelsP[!(levelsP %in% levelsS)]
    })
    hasNewLevels <- sapply(newLevels, length) > 0
    excludeLevels <- any(hasNewLevels)
    time.levs <- Sys.time() - time.levs

    # generate values of new variable
    tt <- Sys.time()
    if ( parallel ) {
      # make cluster        
      cl <- makePSOCKcluster(nr_cores)        
      setDefaultCluster(cl)  
      clusterEvalQ(cl, {"generate.values"})
      clusterExport(cl=cl, varlist=c("generate.values","dataS","strata",
        "dataP","predNames","indStrata","excludeLevels","hasNewLevels",
        "newLevels","formula.cmd","multinom","method","eps","i",
        "limit","censor","spSample","levelsResponse"))        

      values <- parLapply(cl, levels(dataS[, strata]), function(x) {
        generate.values(
          dataSample=dataS[dataS[, strata] == x, , drop=FALSE], 
          dataPop=dataP[indStrata[[x]], predNames, drop=FALSE],
          method=method,
          cur.var=i,
          excludeLevels=excludeLevels, 
          hasNewLevels=hasNewLevels, 
          newLevels=newLevels,
          w=w,
          formula.cmd=formula.cmd,
          eps=eps,
          limit=limit,
          censor=censor,
          levelsResponse=levelsResponse
        )          
      })
      stopCluster(cl)        
    } else {
      values <- lapply(levels(dataS[, strata]), function(x) { 
        generate.values(
          dataSample=dataS[dataS[, strata] == x, , drop=FALSE], 
          dataPop=dataP[indStrata[[x]], predNames, drop=FALSE],
          method=method,
          cur.var=i,
          excludeLevels=excludeLevels, 
          hasNewLevels=hasNewLevels, 
          newLevels=newLevels,
          w=w,
          formula.cmd=formula.cmd,
          eps=eps,
          limit=limit,
          censor=censor,
          levelsResponse=levelsResponse
        )
      })  
    }

    tt <- Sys.time() - tt
    cat("time.levs: "); print(time.levs)
    cat("time.values: "); print(tt)

    values <- factor(unsplit(values, dataP[, strata, drop=FALSE]), levels=levelsResponse)

    ## add new categorical variable to data set
    dataP[, i] <- values
    predNames <- c(predNames, i)
  }

  # return simulated data
  dataP
}


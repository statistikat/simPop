# cur.var: current.variable (value of 'i' from outer loop)
generate.values <- function(dataSample, dataPop, params) {
  if ( !nrow(dataSample) ) {
    return(character())
  }

  meth <- params$method
  cur.var <- params$cur.var
  excludeLevels <- params$excludeLevels
  hasNewLevels <- params$hasNewLevels
  newLevels <- params$newLevels
  w <- params$w
  formula.cmd <- params$formula.cmd
  eps <- params$eps
  limit <- params$limit
  censor <- params$censor
  levelsResponse <- params$levelsResponse  

  indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
  grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]

  # in sample, observations with NAs have been removed to fit the 
  # model, hence population can have additional levels
  # these need to be removed since those probabilities cannot 
  # be predicted from the model
  if ( excludeLevels ) {
    exclude <- mapply(function(pop, new) pop %in% new, 
      pop=grid[, hasNewLevels, drop=FALSE], 
      new=newLevels[hasNewLevels])
    if ( is.null(dim(exclude)) ) {
      exclude <- which(any(exclude))
    } else exclude <- which(apply(exclude, 1, any))
  } else {
    exclude <- integer()
  }
  # fit multinomial model
  # command needs to be constructed as string
  # this is actually a pretty ugly way of fitting the model
  #browser()
  mod <- eval(parse(text=formula.cmd))  # fitted model

  # predict probabilities
  if( length(exclude) == 0 ) {
    newdata <- grid
  } else {
    newdata <- grid[-exclude, , drop=FALSE]
  }
  ind <- match(colnames(newdata), colnames(dataSample))
  for ( i in 1:length(ind) ) {
    if ( is.factor(newdata[,i]) ) {
      newdata[,i] <- factor(as.character(newdata[,i]), levels(dataSample[,ind[i]]))
    }
  }   
  
  if ( meth %in% "multinom" ) {
    probs <- predict(mod, newdata=newdata, type="probs")
  }
  if ( meth %in% "naivebayes" ) {
    probs <- predict(mod, newdata=newdata, type="raw")
  }        
  # TODO: fix error if level-sets are not equal!
  if ( meth %in% "ctree" ) {
    probs <- do.call("rbind", predict(mod, newdata=newdata, type="prob"))
  }
  if ( meth == "liblinear" ) {
    probs <- predict(mod, newx=newdata, proba=TRUE)
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
  if ( length(ind) == 1 ) {
    resample <- function(k, n, p) rep.int(1, n[k])
  } else if ( length(ind) == 2 ) {
    resample <- function(k, n, p) spSample(n[k], c(1-p[k],p[k]))
  } else {
    resample <- function(k, n, p) spSample(n[k], p[k,])
  }
  # generate realizations for each combination
  if ( length(exclude) == 0 ) {
    ncomb <- as.integer(sapply(indGrid, length))
    sim <- lapply(1:length(ncomb), resample, ncomb, probs)
  } else {
    ncomb <- as.integer(sapply(indGrid[-exclude], length))
    sim <- as.list(rep.int(NA, length(indGrid)))
    sim[-exclude] <- lapply(1:length(ncomb), resample, ncomb, probs)
  }
  sim <- unsplit(sim, dataPop, drop=TRUE)
  invisible(levelsResponse[ind][sim])
}

# simulation of variables using random draws from the observed 
# conditional distributions of their multivariate realizations
generate.values.distribution <- function(dataSample, dataPop, params) {
  grid <- params$grid
  additional <- params$additional
  basic <- params$basic
  w <- params$w

  if( !nrow(dataSample) ) {
    return(character())
  }
  
  # population data
  splitS <- split(1:nrow(dataSample), dataSample[, basic], drop=TRUE)
  pSplit <- lapply(splitS, function(i) {
    tmp <- tableWt(dataSample[i, additional], dataSample[i, w])
    tmp <- as.data.frame(tmp)
    p <- ncol(tmp)
    tmp[, p]/sum(tmp[, p])
  })
  splitP <- split(1:nrow(dataPop), dataPop[, basic])
  NSplit <- sapply(splitP, length)
  # in sample, observations with NAs have been removed to fit the 
  # model, hence population can have additional levels
  whichP <- which(names(splitP) %in% names(splitS))
  # generate realizations for each combination
  sim <- as.list(rep.int(NA, length(splitP)))
  sim[whichP] <- mapply(spSample, NSplit[whichP], pSplit, SIMPLIFY=FALSE)
  sim <- unsplit(sim, dataPop[, basic])
  sim <- grid[sim,]
  rownames(sim) <- rownames(dataPop)
  sim
}

simCategorical <- function(dataS, dataP, w = "rb050", strata = "db040",
  basic, additional = c("pl030", "pb220a"),
  method = c("multinom", "distribution","ctree","naivebayes","liblinear"), 
  limit = NULL, censor = NULL, maxit = 500, 
  MaxNWts = 1500, eps = NULL, seed) {
  
  parallel <- FALSE
  if ( Sys.info()["sysname"] != "Windows" ) {
    nr_cores <- detectCores()
    if ( nr_cores > 2 ) {
      parallel <- TRUE
      nr_cores <- nr_cores-1 # keep one core available
    } else {
      parallel <- FALSE
    }  
  }  
  ##### initializations
  if ( !missing(seed) ) {
    set.seed(seed)  # set seed of random number generator
  }
  if ( length(strata) != 1 ) { 
    stop("currently 'strata' must specify exactly one column of 'data'")
  }
  method <- match.arg(method)  
  if ( method == "ctree" ) {
    stop("fix errow with unequal levels of factors in train and test-dataset!\n")
  }  
  
  if ( missing(basic) ) {
    basic <- c("age", "rb090")
    if ( method == "multinom" ) {
      basic <- c(basic, "hsize")
    }
  }
  varNames <- c(w=w, strata=strata, basic, additional)

  # check data
  if ( all(varNames %in% names(dataS)) ) {
    dataS <- dataS[, varNames]
  } else {
    stop("undefined variables in the sample data")
  }
  if ( !all(c(strata, basic) %in% names(dataP)) ) {
    stop("undefined variables in the population data")
  }

  # observations with missings are excluded from simulation
  exclude <- getExclude(dataS)
  if ( length(exclude) ) {
    dataS <- dataS[-exclude,]
  }
  # variables are coerced to factors
  dataS <- checkFactor(dataS, c(strata, basic, additional))
  dataP <- checkFactor(dataP, c(strata, basic))

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

  # list indStrata contains the indices of dataP split by strata
  N <- nrow(dataP)
  indStrata <- split(1:N, dataP[, strata, drop=FALSE])

  ##### simulation
  # predictor variables
  predNames <- basic  # names of predictor variables
  
  if ( method == "distribution" ) {
    params <- list()
    params$grid <- expand.grid(lapply(dataS[, additional], levels))
    params$additional <- additional
    params$basic <- basic
    params$w <- w

    if ( parallel ) {      
      values <- mclapply(levels(dataS[, strata]), function(x) { 
        generate.values.distribution(
          dataSample=dataS[dataS[, strata] == x, , drop=FALSE],
          dataPop=dataP[indStrata[[x]], , drop=FALSE], params)
      })
    } else {
      values <- lapply(levels(dataS[, strata]), function(x) { 
        generate.values.distribution(
          dataSample=dataS[dataS[, strata] == x, , drop=FALSE], 
          dataPop=dataP[indStrata[[x]], , drop=FALSE], params)
      })
    }
    values <- unsplit(values, dataP[, strata, drop=FALSE])

    ## add new categorical variables to data set and return
    dataP[, additional] <- values
    return(dataP)
  }

  # any other method
  for ( i in additional ) {
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
    if ( method == "ctree" ) {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("suppressWarnings(ctree(", formula.cmd, ", weights=as.integer(dataSample$", w, "), data=dataSample))", sep="")						
    }
    if ( method == "naivebayes" ) {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("naiveBayes(", formula.cmd, ", data=dataSample, usekernel=TRUE)", sep="")						
    }
    if ( method == "liblinear" ) {
      formula.cmd <- paste("LiblineaR(data=as.matrix(dataSample[,c('",paste(predNames, collapse="','"),"')]), labels=dataSample[,cur.var], type=0, bias=TRUE,verbose=FALSE)", collapse="", sep="")						
    }

    # check if population data contains factor levels that do not exist 
    # in the sample
    newLevels <- lapply(predNames, function(nam) {
      levelsS <- levels(dataS[, nam])
      levelsP <- levels(dataP[, nam])
      levelsP[!(levelsP %in% levelsS)]
    })
    hasNewLevels <- sapply(newLevels, length) > 0
    excludeLevels <- any(hasNewLevels)

    # generate values of new variable
    params <- list()
    params$method <- method
    params$cur.var <- i
    params$excludeLevels <- excludeLevels
    params$hasNewLevels <- hasNewLevels
    params$newLevels <- newLevels
    params$w <- w
    params$formula.cmd <- formula.cmd
    params$eps <- eps
    params$limit <- limit
    params$censor <- censor
    params$levelsResponse <- levelsResponse
    if ( parallel ) {
      values <- mclapply(levels(dataS[, strata]), function(x) {
        generate.values(
          dataSample=dataS[dataS[, strata] == x, , drop=FALSE],
          dataPop=dataP[indStrata[[x]], predNames, drop=FALSE], params
        )        
      })
    } else {
      values <- lapply(levels(dataS[, strata]), function(x) {
        generate.values(
          dataSample=dataS[dataS[, strata] == x, , drop=FALSE],
          dataPop=dataP[indStrata[[x]], predNames, drop=FALSE], params
        )
      })
    }
    values <- factor(unsplit(values, dataP[, strata, drop=FALSE]), levels=levelsResponse)
    ## add new categorical variable to data set
    dataP[, i] <- values
    predNames <- c(predNames, i)
  }
  invisible(dataP)
}

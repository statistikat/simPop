generateValues <- function(dataSample, dataPop, params) {
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
  grid <- as.data.frame(grid)

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
    } else {
      exclude <- which(apply(exclude, 1, any))
    }
  } else {
    exclude <- integer()
  }
  # fit multinomial model
  # command needs to be constructed as string
  # this is actually a pretty ugly way of fitting the model
  mod <- eval(parse(text=formula.cmd))  # fitted model

  # predict probabilities
  if ( length(exclude) == 0 ) {
    newdata <- grid
  } else {
    newdata <- grid[-exclude, , drop=FALSE]
  }
  ind <- match(colnames(newdata), colnames(dataSample))
  for ( i in 1:length(ind) ) {
    if ( is.factor(newdata[,i]) ) {
      newdata[,i] <- factor(as.character(newdata[,i]), levels(dataSample[[ind[i]]]))
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
  # set too small probabilities to exactly 0
  if ( !is.null(eps) ) {
    probs[probs < eps] <- 0
  }

  # ensure code works for missing levels of response
  ind <- as.integer(which(table(dataSample[[cur.var]]) > 0))
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
generateValues_distribution <- function(dataSample, dataPop, params) {
  grid <- params$grid
  additional <- params$additional
  basic <- params$basic
  w <- params$w

  if( !nrow(dataSample) ) {
    return(character())
  }

  # population data
  splitS <- split(1:nrow(dataSample), dataSample[, basic, with=F], drop=TRUE)
  pSplit <- lapply(splitS, function(i) {
    tmp <- tableWt(dataSample[i, additional, with=F], dataSample[[w]][i])
    tmp <- as.data.frame(tmp)
    p <- ncol(tmp)
    tmp[, p]/sum(tmp[, p])
  })
  splitP <- split(1:nrow(dataPop), dataPop[, basic, with=F])
  NSplit <- sapply(splitP, length)
  # in sample, observations with NAs have been removed to fit the
  # model, hence population can have additional levels
  whichP <- which(names(splitP) %in% names(splitS))
  # generate realizations for each combination
  sim <- as.list(rep.int(NA, length(splitP)))
  sim[whichP] <- mapply(spSample, NSplit[whichP], pSplit, SIMPLIFY=FALSE)
  sim <- unsplit(sim, dataPop[, basic, with=F])
  sim <- grid[sim,]
  rownames(sim) <- rownames(dataPop)
  sim
}

simCategorical <- function(synthPopObj, additional,
  method=c("multinom", "distribution", "naivebayes"),
  limit=NULL, censor=NULL, maxit=500, MaxNWts=1500, eps=NULL, seed=1) {

  dataP <- synthPopObj@pop
  dataS <- synthPopObj@sample
  data_pop <- dataP@data
  data_sample <- dataS@data
  basic <- synthPopObj@basicHHvars

  if ( any(additional %in% colnames(data_pop)) ) {
    stop("variables already exist in the population!\n")
  }

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

  # observations with missings are excluded from simulation
  exclude <- getExclude(data_sample)
  if ( length(exclude) > 0 ) {
    data_sample <- data_sample[-exclude,]
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

  # list indStrata contains the indices of dataP split by strata
  N <- nrow(data_pop)
  indStrata <- split(1:N, data_pop[[dataP@strata]])

  ##### simulation
  # predictor variables
  predNames <- basic  # names of predictor variables

  if ( method == "distribution" ) {
    params <- list()
    params$grid <- expand.grid(lapply(data_sample[,additional, with=F], levels))
    params$additional <- additional
    params$basic <- basic
    params$w <- dataS@weight

    if ( parallel ) {
      values <- mclapply(levels(data_sample[[dataS@strata]]), function(x) {
        generateValues_distribution(
          dataSample=data_sample[data_sample[[dataS@strata]] == x,],
          dataPop=data_pop[indStrata[[x]], predNames, with=F], params)
      })
    } else {
      values <- lapply(levels(data_sample[[dataS@strata]]), function(x) {
        generateValues_distribution(
          dataSample=data_sample[data_sample[[dataS@strata]] == x,],
          dataPop=data_pop[indStrata[[x]], predNames, with=F], params)
      })
    }
    values <- do.call("rbind", values)

    ## add new categorical variables to data set and return
    for ( i in additional ) {
      data_pop[[i]] <- values[,i]
    }
    synthPopObj@pop@data <- data_pop
    return(invisible(synthPopObj))
  }

  # any other method
  for ( i in additional ) {
    # components of multinomial model are specified
    levelsResponse <- levels(data_sample[[i]])

    # simulation of variables using a sequence of multinomial models
    if ( method == "multinom" ) {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("suppressWarnings(multinom(", formula.cmd,
        ", weights=", dataS@weight, ", data=dataSample, trace=FALSE",
        ", maxit=",maxit, ", MaxNWts=", MaxNWts,"))", sep="")
    }
    # simulation via recursive partitioning and regression trees
    if ( method == "ctree" ) {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("suppressWarnings(ctree(", formula.cmd, ", weights=as.integer(dataSample$", dataS@weight, "), data=dataSample))", sep="")
    }
    if ( method == "naivebayes" ) {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd <- paste("naiveBayes(", formula.cmd, ", data=dataSample, usekernel=TRUE)", sep="")
    }

    # check if population data contains factor levels that do not exist
    # in the sample
    newLevels <- lapply(predNames, function(nam) {
      levelsS <- levels(data_sample[[nam]])
      levelsP <- levels(data_pop[[nam]])
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
    params$w <- dataS@weight
    params$formula.cmd <- formula.cmd
    params$eps <- eps
    params$limit <- limit
    params$censor <- censor
    params$levelsResponse <- levelsResponse
    if ( parallel ) {
      values <- mclapply(levels(data_sample[[dataS@strata]]), function(x) {
        generateValues(
          dataSample=data_sample[data_sample[[dataS@strata]] == x,],
          dataPop=data_pop[indStrata[[x]], predNames, with=F], params
        )
      })
    } else {
      values <- lapply(levels(data_sample[[dataS@strata]]), function(x) {
        generateValues(
          dataSample=data_sample[data_sample[[dataS@strata]] == x,],
          dataPop=data_pop[indStrata[[x]], predNames, with=F], params
        )
      })
    }
    values <- factor(unsplit(values, data_pop[[dataP@strata]]), levels=levelsResponse)
    ## add new categorical variable to data set
    data_pop[[i]] <- values
    predNames <- c(predNames, i)
  }
  synthPopObj@pop@data <- data_pop
  invisible(synthPopObj)
}

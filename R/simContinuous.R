generateValues_multinom <- function(dataSample, dataPop, params) {
  excludeLevels <- params$excludeLevels
  maxit <- params$maxit
  MaxNWts <- params$MaxNWts
  command <- params$command
  eps <- params$eps
  limit <- params$limit
  censor <- params$censor
  hasNewLevels  <- params$hasNewLevels
  newLevels <- params$newLevels
  name <- params$name
  response <- params$response

  # unique combinations in the stratum of the population need
  # to be computed for prediction
  indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
  grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]
  grid <- as.data.frame(grid)

  # in sample, observations with NAs have been removed to fit the
  # model, hence population can have additional levels
  # these need to be removed since those probabilities cannot
  # be predicted from the model
  if ( excludeLevels ) {
    exclude <- mapply(function(pop, new) pop %in% new,
      pop=grid[, hasNewLevels, drop=FALSE], new=newLevels[hasNewLevels]
    )
    if ( is.null(dim(exclude)) ) {
      exclude <- which(any(exclude))
    } else {
      exclude <- which(apply(exclude, 1, any))
    }
  } else {
    exclude <- integer()
  }
  # fit multinomial model
  mod <- eval(parse(text=command))  # fitted model
  # predict probabilities
  if ( length(exclude) == 0 ) {
    probs <- predict(mod, newdata=grid, type="probs")
  } else {
    probs <- predict(mod, newdata=grid[-exclude, , drop=FALSE], type="probs")
  }
  # set too small probabilities to exactly 0
  if ( !is.null(eps) ) {
    probs[probs < eps] <- 0
  }
  # ensure it works for missing levels of response
  ind <- as.integer(which(table(dataSample[[name]]) > 0))
  if ( length(ind) > 2 && (nrow(grid)-length(exclude)) == 1 ) {
    probs <- t(probs)
  }
  # account for structural zeros
  if ( (!is.null(limit) || !is.null(censor)) && !is.null(dim(probs)) ) {
    if ( length(exclude) == 0 ) {
      probs <- adjustProbs(probs, grid, names(indGrid), limit, censor)
    } else {
      probs <- adjustProbs(probs, grid[-exclude, , drop=FALSE], names(indGrid)[-exclude], limit, censor)
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
  sim <- as.list(rep.int(NA, length(indGrid)))
  if ( length(exclude) == 0 ) {
    ncomb <- as.integer(sapply(indGrid, length))
    sim <- lapply(1:length(ncomb), resample, ncomb, probs)
  } else {
    ncomb <- as.integer(sapply(indGrid[-exclude], length))
    sim[-exclude] <- lapply(1:length(ncomb), resample, ncomb, probs)
  }
  sim <- unsplit(sim, dataPop, drop=TRUE)
  # return realizations
  levels(response)[ind][sim]
}

generateValues_lm <- function(dataSample, dataPop, params) {
  if ( !nrow(dataSample) ) {
    return(numeric())
  }
  coef <- params$coef
  excludeLevels <- params$excludeLevels
  hasNewLevels <- params$hasNewLevels
  newLevels <- params$newLevels
  command <- params$command
  predNames <- params$predNames
  additional <- params$additional
  const <- params$const
  formula <- params$formula
  levels <- params$levels
  residuals <- params$residuals
  log <- params$log

  # unique combinations in the stratum of the population need to be computed for prediction
  indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
  grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]
  grid <- as.data.frame(grid)

  # in sample, observations with NAs have been removed to fit the
  # model, hence population can have additional levels
  # these need to be removed since those probabilities cannot
  # be predicted from the model
  if ( excludeLevels ) {
    exclude <- mapply(function(pop, new) pop %in% new,
      pop=grid[, hasNewLevels, drop=FALSE], new=newLevels[hasNewLevels]
    )
    if ( is.null(dim(exclude)) ) {
      exclude <- which(any(exclude))
    } else {
      exclude <- which(apply(exclude, 1, any))
    }
    if ( length(exclude) > 0 ) {
      grid <- grid[-exclude, , drop=FALSE]
    }
    for ( j in predNames[hasNewLevels] ) {
      # drop new factor levels
      grid[, j] <- factor(as.character(grid[, j]), levels=levels(dataSample[[j]]))
    }
  }
  # fit linear model
  mod <- eval(parse(text=command))
  # add coefficients from auxiliary model if necessary
  tmp <- coef
  coef[names(coef(mod))] <- coef(mod)
  mod$coefficients <- coef
  # prediction
  # add 0 variable to combinations for use of 'model.matrix'
  newdata <- cbind(grid, 0)
  names(newdata) <- c(predNames, additional)
  newdata <- model.matrix(formula, data=newdata)
  if ( length(exclude) == 0 ) {
    pred <- spPredict(mod, newdata)
  } else {
    pred <- as.list(rep.int(NA, length(indGrid)))
    pred[-exclude] <- spPredict(mod, newdata)
  }
  pred <- unsplit(pred, dataPop, drop=TRUE)
  # add error terms
  if ( residuals ) {
    error <- sample(residuals(mod), size=nrow(dataPop), replace=TRUE)
  } else {
    mu <- median(residuals(mod))
    sigma <- mad(residuals(mod))
    error <- rnorm(nrow(dataPop), mean=mu, sd=sigma)
  }
  # return realizations
  sim <- pred + error
  if ( log ) {
    res <- exp(sim)  # transform back
    if ( !is.null(const) ) {
      res <- res - const  # subtract constant
    }
    return(res)
  } else {
    return(sim)
  }
}

generateValues_binary <- function(dataSample, dataPop, params) {
  excludeLevels <- params$excludeLevels
  hasNewLevels <- params$hasNewLevels
  newLevels <- params$newLevels
  predNames <- params$predNames
  name <- params$name
  weight <- params$weight
  useAux <- params$useAux
  tol <- params$tol
  eps <- params$eps
  if ( !nrow(dataSample) ) {
    return(numeric())
  }

  # unique combinations in the stratum of the population need to be computed for prediction
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
      new=newLevels[hasNewLevels]
    )
    if ( is.null(dim(exclude)) ) {
      exclude <- which(any(exclude))
    } else {
      exclude <- which(apply(exclude, 1, any))
    }
    if ( length(exclude) > 0 ) {
      grid <- grid[exclude, , drop=FALSE]
    }
    for ( j in predNames[hasNewLevels] ) {
      # drop new factor levels
      grid[, j] <- factor(as.character(grid[, j]), levels=levels(dataSample[[j]]))
    }
  }
  # add 0 variable to combinations for use of 'model.matrix'
  Xnew <- cbind(grid, 0)
  names(Xnew) <- c(predNames, name)
  Xnew <- model.matrix(params$command, data=Xnew)

  # fit logit model
  X <- model.matrix(params$command, data=dataSample)
  y <- dataSample[[name]]
  weights <- dataSample[[weight]]
  mod <- logitreg(X, y, weights=weights)
  # add parameters from auxiliary model if necessary
  if ( useAux ) {
    indPar <- abs(mod$par) < tol
    mod$par[indPar] <- params$par[indPar]
  }
  # predict probabilities
  tmp <- exp(Xnew %*% mod$par)
  # avoid integer overflow
  p <- ifelse(is.infinite(tmp), 1, as.numeric(tmp / (1 + tmp)))
  # set too small probabilities to exactly 0
  if ( !is.null(eps) ) {
    p[p < eps] <- 0
  }
  # generate realizations for each combination
  if ( length(exclude) == 0 ) {
    ncomb <- as.integer(sapply(indGrid, length))
    sim <- lapply(1:length(ncomb), function(k) {
      spSample(ncomb[k], c(1-p[k], p[k])) - 1
    })
  } else {
    ncomb <- as.integer(sapply(indGrid[-exclude], length))
    sim <- as.list(rep.int(NA, length(indGrid)))
    sim[-exclude] <- lapply(1:length(ncomb), function(k) {
      spSample(ncomb[k], c(1-p[k], p[k])) - 1
    })
  }
  # return realizations
  unsplit(sim, dataPop, drop=TRUE)
}

runModel <- function(dataS, dataP, params, typ) {
  x <- NULL
  strata <- params$strata
  pp <- parallelParameters(nr_cpus=params$nr_cpus, nr_strata=length(levels(dataS[[strata]])))
  indStrata <- params$indStrata
  predNames <- params$predNames
  additional <- c(params$additional, params$name)

  if ( pp$parallel ) {
    # windows
    if ( pp$have_win ) {
      cl <- makePSOCKcluster(pp$nr_cores)
      registerDoParallel(cl,cores=pp$nr_cores)
      if ( typ=="multinom" ) {
        valuesCat <- foreach(x=levels(dataS[[strata]]), .options.snow=list(preschedule=TRUE)) %dopar% {
          generateValues_multinom(
            dataSample=dataS[dataS[[strata]] == x,],
            dataPop=dataP[indStrata[[x]], predNames, with=F], params
          )
        }
      }
      if ( typ=="binary" ) {
        valuesCat <- foreach(x=levels(dataS[[strata]]), .options.snow=list(preschedule=TRUE)) %dopar% {
          generateValues_binary(
            dataSample=dataS[dataS[[strata]] == x,],
            dataPop=dataP[indStrata[[x]], predNames, with=F], params
          )
        }
      }
      if ( typ=="lm" ) {
        valuesCat <- foreach(x=levels(dataS[[strata]]), .options.snow=list(preschedule=TRUE)) %dopar% {
          generateValues_lm(
            dataSample=dataS[dataS[[strata]] == x,],
            dataPop=dataP[indStrata[[x]], predNames, with=F], params
          )
        }
      }
      stopCluster(cl)
    }

    # linux/mac
    if ( !pp$have_win ) {
      if ( typ=="multinom" ) {
        valuesCat <- mclapply(levels(dataS[[strata]]), function(x) {
          generateValues_multinom(
            dataSample=dataS[dataS[[strata]] == x,],
            dataPop=dataP[indStrata[[x]], predNames, with=F], params
          )
        },mc.cores=pp$nr_cores)
      }
      if ( typ=="binary" ) {
        valuesCat <- mclapply(levels(dataS[[strata]]), function(x) {
          generateValues_binary(
            dataSample=dataS[dataS[[strata]] == x,],
            dataPop=dataP[indStrata[[x]], predNames, with=F], params
          )
        },mc.cores=pp$nr_cores)
      }
      if ( typ=="lm" ) {
        valuesCat <- mclapply(levels(dataS[[strata]]), function(x) {
          generateValues_lm(
            dataSample=dataS[dataS[[strata]] == x,],
            dataPop=dataP[indStrata[[x]], predNames, with=F], params
          )
        },mc.cores=pp$nr_cores)
      }
    }
  } else {
    if ( typ=="multinom" ) {
      valuesCat <- lapply(levels(dataS[[strata]]), function(x) {
        generateValues_multinom(
          dataSample=dataS[dataS[[strata]] == x,c(predNames, additional), with=F],
          dataPop=dataP[indStrata[[x]], predNames, with=F], params
        )
      })
    }
    if ( typ=="binary" ) {
      valuesCat <- lapply(levels(dataS[[strata]]), function(x) {
        generateValues_binary(
          dataSample=dataS[dataS[[strata]] == x,],
          dataPop=dataP[indStrata[[x]], predNames, with=F], params
        )
      })
    }
    if ( typ=="lm" ) {
      valuesCat <- lapply(levels(dataS[[strata]]), function(x) {
        generateValues_lm(
          dataSample=dataS[dataS[[strata]] == x,],
          dataPop=dataP[indStrata[[x]], predNames, with=F], params
        )
      })
    }
  }

  # check for errors
  res <- sapply(valuesCat, class)
  if ( any(res=="try-error") ) {
    stop(paste0("Error in estimating the linear model. Try to specify a more simple model!\n"))
  }

  if ( typ=="multinom" ) {
    response <- dataS[[params$name]]
    valuesCat <- factor(unsplit(valuesCat, dataP[[strata]]), levels=levels(response))
  }
  if ( typ=="binary" ) {
    valuesCat <- unsplit(valuesCat, dataP[, strata, drop=FALSE])
  }
  if ( typ=="lm" ) {
    valuesCat <- unlist(valuesCat, dataP[[strata]])
  }
  return(valuesCat)
}

simContinuous <- function(simPopObj, additional = "netIncome",
  method = c("multinom", "lm"), zeros = TRUE,
  breaks = NULL, lower = NULL, upper = NULL,
  equidist = TRUE, probs = NULL, gpd = TRUE,
  threshold = NULL, est = "moments", limit = NULL,
  censor = NULL, log = TRUE, const = NULL,
  alpha = 0.01, residuals = TRUE, keep = TRUE,
  maxit = 500, MaxNWts = 1500,
  tol = .Machine$double.eps^0.5,
  nr_cpus=NULL, eps = NULL, regModel="basic", byHousehold=FALSE, imputeMissings=FALSE, seed) {

  x <- hhid <- vals <- id <- V1 <- NULL

  samp <- simPopObj@sample
  pop <- simPopObj@pop
  basic <- simPopObj@basicHHvars
  strata <- samp@strata
  weight <- samp@weight

  dataS <- samp@data
  dataP <- pop@data

  if ( additional %in% names(dataP)) {
    stop(paste0("Variable '",additional,"' already available in the synthetic population!\n"))
  }

  ## initializations
  if ( !missing(seed) ) {
    set.seed(seed)
  }

  if ( length(additional) != 1 ) {
    stop("currently only one additional variable can be generated at a time")
  }
  if ( !additional %in% colnames(samp@data) ) {
    stop("variable 'additional' must be included in the sample of input 'simPopObj'!\n")
  }

  regInput <- regressionInput(simPopObj, additional=additional, regModel=regModel)
  predNames <- regInput[[1]]$predNames
  estimationModel <- regInput[[1]]$formula

  varNames <- unique(c(predNames, weight, additional, strata))
  dataS <- dataS[,varNames, with=F]

  method <- match.arg(method)
  zeros <- isTRUE(zeros)
  log <- isTRUE(log)

  if ( is.numeric(alpha) && length(alpha) > 0 ) {
    alpha <- rep(alpha, length.out=2)
    if ( !all(is.finite(alpha)) || any(alpha < 0) || sum(alpha) >= 1 ) {
      alpha <- NULL
      warning("invalid parameter 'alpha': trimming is not applied\n")
    }
  } else {
    alpha <- NULL
  }

  # observations with missings are excluded from simulation
  #exclude <- getExclude(dataS[,c(additional,predNames),with=F]) # fixes #31?
  #if ( length(exclude) ) {
  #  dataS <- dataS[-exclude,]
  #}

  # temporarily impute (using hotdeck) / or check (if imputeMissings=FALSE)
  # missing values in additional variables in the sample
  if ( is.null(samp@strata) ) {
    modelVars <- setdiff(predNames, c(weight,basic,pop@hhsize))
  } else {
    modelVars <- setdiff(predNames, c(strata,weight,basic,pop@hhsize))
  }

  if ( length(modelVars) > 0 & imputeMissings ) {
    dataS_orig <- dataS[,modelVars,with=F]
    dataS <- hotdeck(dataS, variable=modelVars, domain_var=sample@strata, imp_var=FALSE)
  }

  # check for NAs and warn user
  if ( !imputeMissings) {
    naTab <- dataS[,lapply(.SD, is.na), .SDcols=c(additional,predNames)]
    perc.miss <- sum(rowSums(naTab)!=0) / nrow(dataS)
    if ( perc.miss > 0 ) {
      wm <- paste0("There are ~",formatC(100*perc.miss,format="f", digits=1),"% ")
      wm <- paste0(wm, "observations in the response/predictors with at least one missing variable.\n")
      wm <- paste0(wm, "If you get errors in the estimation procedure, consider to recode these missing ")
      wm <- paste0(wm, "values (e.g. by assigning an additional category) or try to specify a different model.\n\n")
      for ( z in 1:ncol(naTab) ) {
        vv <- colnames(naTab)[z]
        missv <- sum(naTab[[z]])
        missp <- formatC(100*missv/nrow(dataS),format="f", digits=1)
        if ( vv == additional ) {
          wm <- paste0(wm, "Variable '",vv,"' (response): ",missv," missing values (~",missp,"%).\n")
        } else {
          wm <- paste0(wm, "Variable '",vv,"' (predictor): ",missv," missing values (~",missp,"%).\n")
        }
      }
      warning(wm)
    }
  }

  # variables are coerced to factors
  select <- unique(c(predNames, samp@strata)) # strata always included
  dataS <- checkFactor(dataS, select)
  dataP <- checkFactor(dataP, select)

  # sample data of variable to be simulated
  additionalS <- dataS[[additional]]

  ## determine which models to fit and do further initializations
  haveBreaks <- !is.null(breaks)
  if ( method == "multinom" ) {
    useMultinom <- TRUE
    useLogit <- FALSE
    useLm <- FALSE
    # define break points (if missing)
    if ( haveBreaks ) {
      checkBreaks(breaks)
      breaks <- if(zeros) union(breaks, 0) else unique(breaks)
      breaks <- sort(breaks)
    } else {
      if ( is.null(upper) && gpd ) {
        upper <- Inf
      }
      breaks <- getBreaks(additionalS, dataS[[weight]], zeros, lower, upper, equidist, probs)
    }
  } else {
    useLm <- TRUE
    if ( log ) {
      if ( is.null(const) ) {
        ## use log-transformation
        # check for negative values
        neg <- which(additionalS < 0)
        haveNeg <- length(neg) > 0
        if ( haveNeg ) {
          # define break points for negative values
          if ( haveBreaks ) {
            checkBreaks(breaks)
            breaks <- c(unique(breaks[breaks < 0]), 0)
          } else {
            breaks <- getBreaks(additionalS[neg], dataS[[weight]][neg], zeros=TRUE, lower, upper)
          }
          if ( zeros || length(breaks) > 2 ) {
            useMultinom <- TRUE
            breaks <- c(breaks, Inf)  # add Inf to breakpoints
          } else {
            useMultinom <- FALSE
          }
          useLogit <- !useMultinom
        } else {
          useLogit <- zeros || any(additionalS == 0)
          useMultinom <- FALSE
        }
      } else {
        # check constant
        if ( !is.numeric(const) || length(const) == 0 ) {
          stop("'const' must be numeric\n")
        } else {
          const <- const[1]
        }
        # set control parameters
        useLogit <- zeros || any(additionalS == 0)
        useMultinom <- FALSE
      }
    } else {
      # logistic model is used in case of semi-continuous variable
      useLogit <- zeros
      # multinomial model is not needed
      useMultinom <- FALSE
    }
  }

  ## some general preparations for the simulation
  # list indStrata contains the indices of dataP split by strata
  N <- nrow(dataP)
  indP <- 1:N
  indStrata <- split(indP, dataP[[strata]])
  #fpred <- paste(predNames, collapse = " + ")  # for formula
  # check if population data contains factor levels that do not exist
  # in the sample
  newLevels <- lapply(predNames, function(nam) {
    levelsS <- levels(dataS[[nam]])
    levelsP <- levels(dataP[[nam]])
    levelsP[!(levelsP %in% levelsS)]
  })
  hasNewLevels <- sapply(newLevels, length) > 0
  excludeLevels <- any(hasNewLevels)

  ## preparations for multinomial or binomial logit model
  if ( useMultinom || useLogit ) {
    name <- getCatName(additional)
    estimationModel <- gsub(additional, name, estimationModel)
  }

  if ( useMultinom ) {
    ## some preparations
    dataS[[name]] <- getCat(additionalS, breaks, zeros, right=TRUE)
    response <- dataS[[name]]  # response variable
    # check threshold for GPD (if supplied)
    if ( !useLm && gpd && !is.null(threshold) && length(threshold) != 1 ) {
      stop("'threshold' must be a single numeric value")
    }

    ## simulate categories
    # TODO: share code with 'simCategorical'
    params <- list()
    params$excludeLevels <- excludeLevels
    # command needs to be constructed as string
    # this is actually a pretty ugly way of fitting the model
    params$command <- paste("suppressWarnings(multinom(", estimationModel,
      ", weights=", weight, ", data=dataSample, trace=FALSE",
      ", maxit=maxit, MaxNWts=MaxNWts))", sep="")
    params$maxit <- maxit
    params$MaxNWts <- MaxNWts
    params$eps <- eps
    params$limit <- limit
    params$censor <- censor
    params$hasNewLevels  <- hasNewLevels
    params$newLevels <- newLevels
    params$name <- name
    params$response <- response
    params$strata <- strata
    params$nr_cpus <- nr_cpus
    params$indStrata <- indStrata
    params$predNames <- predNames
    params$additional <- c(additional, weight)
    cat("running multinom with the following model:\n")
    cat(gsub("))",")",gsub("suppressWarnings[(]","",params$command)),"\n")

    # run in parallel if possible
    valuesCat <- runModel(dataS, dataP, params, typ="multinom")

    ## simulate (semi-)continuous values
    tcat <- table(valuesCat)
    ncat <- length(tcat)

    icat <- 1:ncat
    values <- as.list(rep.int(NA, ncat))
    # zeros

    if ( zeros ) {
      # bug: missing 0 even though zeros is not null?
      izero <- which(breaks == 0)
      values[izero] <- 0
      tcat <- tcat[-izero]
      ncat <- length(tcat)
      icat <- icat[-izero]
    }

    # values to be simulated with linear model or draws from Pareto
    # distribution
    if ( useLm ) {
      # last breakpoint is Inf, the one before is 0
      nunif <- ncat - 1  # leave category of positive values out
    } else {
      nbreaks <- length(breaks)
      if ( gpd ) {
        if ( is.null(threshold) ) {
          if ( !haveBreaks && (!isTRUE(equidist) || !is.null(probs)) ) {
            ngpd <- nbreaks-2
          } else {
            ngpd <- nbreaks-1
          }
        } else if ( any(tmp <- breaks >= threshold) ) {
          ngpd <- min(which(tmp))
        } else {
          ngpd <- nbreaks
        }
      } else {
        ngpd <- nbreaks
      }
      if ( gpd && ngpd <= ncat ) {
        # adjust threshold and fit GPD
        threshold <- breaks[ngpd]  # adjust threshold
        estPar <- fitgpd(additionalS, threshold, est)  # fit GPD
        estPar <- estPar[["fitted.values"]]  # parameters of GPD
        # generalized pareto distribution
        igpd <- ngpd:ncat
        values[icat[igpd]] <- lapply(igpd, function(i) {
          truncPareto(tcat[i], loc=threshold, scale=estPar["scale"], shape=estPar["shape"], breaks[i], breaks[i+1])
        })
      }
      nunif <- ngpd - 1
    }
    # uniform distribution
    if ( nunif > 0 ) {
      iunif <- 1:nunif
      values[icat[iunif]] <- lapply(iunif, function(i) {
        runif(tcat[i], breaks[i], breaks[i+1])
      })
    }
    # turn list into vector of values
    values <- unsplit(values, valuesCat)
  }

  if ( useLogit ) {
    ## some preparations
    if ( log && is.null(const) && haveNeg ) {
      indS <- additionalS > 0
    } else {
      indS <- additionalS != 0
    }
    dataS[[name]] <- as.integer(indS)
    estimationModel <- as.formula(estimationModel)  # formula for model
    # auxiliary model for all strata (used in case of empty combinations)
    useAux <- !is.null(tol)
    if ( useAux ) {
      if ( length(tol) != 1 || tol <= 0 ) {
        stop("'tol' must be a single small positive value!\n")
      }
      X <- model.matrix(estimationModel, data=dataS)
      y <- dataS[[name]]
      weights <- dataS[[weight]]
      mod <- logitreg(X, y, weights=weights)
      par <- mod$par
    } else {
      par <- NULL
    }

    ## simulate binary vector
    params <- list()
    params$excludeLevels <- excludeLevels
    params$hasNewLevels <- hasNewLevels
    params$newLevels <- newLevels
    params$predNames <- predNames
    params$tol <- tol
    params$eps <- eps
    params$weight <- weight
    params$useAux <- useAux
    params$name <- name
    params$strata <- strata
    params$nr_cpus <- nr_cpus
    params$indStrata <- indStrata
    params$predNames <- predNames
    params$additional <- additional
    params$par <- par
    params$command <- estimationModel
    # run in parallel if possible
    valuesCat <- runModel(dataS, dataP, params, typ="binary")
  }

  if ( useLm ) {
    ## some preparations
    if ( useMultinom ) {
      catLm <- names(tcat)[ncat]  # category for positive values
      dataS <- dataS[response == catLm, , drop=FALSE]
      indP <- valuesCat == catLm
    } else if( useLogit ) {
      dataS <- dataS[indS, , drop=FALSE]  # select only non-zeros
      indP <- valuesCat == 1  # indicates non-zeros in population
    }
    if ( useMultinom || useLogit ) {
      # adjust population data
      ii <- which(indP)
      if ( length(ii) > 0 ) {
        dataPop <- dataP[which(indP), , drop=FALSE]
      } else {
        dataPop <- dataP
      }
      # list indStrata is adjusted so that it only contains
      # indices of persons in population with non-zero value
      indStrata <- split(1:nrow(dataPop), dataPop[[strata]])
    } else {
      dataPop <- dataP
    }

    ## trim data (if specified)
    if ( !is.null(alpha) ) {
      additionalS <- dataS[[additional]]
      p <- c(alpha[1], 1-alpha[2])
      bounds <- quantileWt(additionalS, dataS[[weight]], p)
      select <- additionalS > bounds[1] & additionalS < bounds[2]
      dataSample <- dataS[select, , drop=FALSE]
      # check if all relevant levels of predictor variables are still
      # contained in sample after trimming
      # if not, trimming is not applied and a warning message is generated
      check <- unlist(sapply(predNames, function(i) {
        table(dataS[[i]]) > 0 & table(dataSample[[i]]) == 0
      }))
      if ( any(check) ) {
        dataSample <- dataS
        warning("trimming could not be applied\n")
      }
    } else {
      dataSample <- dataS
    }

    ## fit linear model
    # formula for linear model
    if ( log ) {
      fname <- paste("log(", additional, if(!is.null(const)) " + const", ")", sep = "")
    } else {
      fname <- additional
    }
    fstring <- paste0(fname, " ~ ", unlist(strsplit(estimationModel,"~"))[2])
    formula <- as.formula(fstring)
    # auxiliary model for all strata (used in case of empty combinations)
    weights <- dataSample[[weight]]
    mod <- lm(formula, weights=weights, data=dataSample)
    coef <- coef(mod)

    # simulate values
    params <- list()
    params$coef <- coef
    params$command <- paste("lm(", fstring,", weights=", weight, ", data=dataSample)", sep="")
    params$name <- fname
    params$excludeLevels <- excludeLevels
    params$hasNewLevels <- hasNewLevels
    params$newLevels <- newLevels
    params$predNames <- predNames
    params$additional <- c(additional, weight)
    params$const <- const
    params$formula <- formula
    params$residuals <- residuals
    params$log <- log
    params$strata <- strata
    params$nr_cpus <- nr_cpus
    params$indStrata <- indStrata
    params$predNames <- predNames
    params$additional <- c(additional, weight)

    valuesTmp <- runModel(dataS, dataP, params, typ="lm")

    ## put simulated values together
    if ( useMultinom ) {
      values[which(indP == 1)] <- valuesTmp
    } else {
      if ( useLogit ) {
        if ( log && is.null(const) && haveNeg ) {
          # only one category for non-positive values (two breakpoints, one of them is 0)
          values <- rep.int(NA, N)
          nonpos <- which(indP == 0)
          values[nonpos] <- runif(length(nonpos), breaks[1], breaks[2])
        } else {
          values <- ifelse(is.na(indP), NA, 0) # only zeros
        }
        values[which(indP == 1)] <- valuesTmp
      } else {
        values <- valuesTmp
      }
    }
  }

  # reset imputed variables in sample
  if ( imputeMissings ) {
    for ( i in 1:ncol(dataS_orig)) {
      cmd <- paste0("dataS[,",colnames(dataS_orig)[i],":=dataS_orig$",colnames(dataS_orig)[i],"]")
      eval(parse(text=cmd))
    }
  }

  # attach new variable(s) to population data
  if ( useMultinom && keep ) {
    dataP[[name]] <- valuesCat
  }

  # calculate mean of new variable by household
  if ( byHousehold ) {
    xx <- data.table(id=1:length(values), hhid=dataP[[pop@hhid]], vals=values)
    setkey(xx, hhid)
    yy <- xx[,mean(vals, na.rm=TRUE), by=key(xx)]
    xx <- merge(xx, yy, all.x=TRUE)
    setkey(xx, id)
    xx[is.nan(V1), V1:=NA]
    values <- xx$V1
  }

  # return simulated data
  dataP[[additional]] <- values
  simPopObj@pop@data <- dataP
  invisible(simPopObj)
}

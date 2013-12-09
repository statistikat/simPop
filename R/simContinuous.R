# ----------------------------------------
# Authors: Andreas Alfons and Stefan Kraft
#          Vienna University of Technology
# ----------------------------------------

# TODO: allow for continuous predictors

simContinuous <- function(dataS, dataP, w = "rb050", strata = "db040",
        basic = c("age", "rb090", "hsize"), 
        additional = "netIncome", 
        method = c("multinom", "lm"), zeros = TRUE, 
        breaks = NULL, lower = NULL, upper = NULL, 
        equidist = TRUE, probs = NULL, gpd = TRUE, 
        threshold = NULL, est = "moments", limit = NULL, 
        censor = NULL, log = TRUE, const = NULL, 
        alpha = 0.01, residuals = TRUE, keep = TRUE, 
        maxit = 500, MaxNWts = 1500, 
        tol = .Machine$double.eps^0.5, 
        eps = NULL, seed) {
    
    ## initializations
    if(!missing(seed)) set.seed(seed)  # set seed of random number generator
    if(length(strata) != 1) {
        stop("currently 'strata' must specify exactly one column of 'data'")
    }
    if(length(additional) != 1) {
        stop("currently only one additional", 
            " variable can be generated at a time")
    }
    varNames <- c(w=w, strata=strata, basic, additional)
    method <- match.arg(method)
    zeros <- isTRUE(zeros)
    log <- isTRUE(log)
    # check data
    if(all(varNames %in% colnames(dataS))) {
		dataS <- dataS[, varNames]
	} else stop("undefined variables in the sample data")
    if(!all(c(strata, basic) %in% colnames(dataP))) {
        stop("undefined variables in the population data")
    }
    if(is.numeric(alpha) && length(alpha) > 0) {
        alpha <- rep(alpha, length.out=2)
        if(!all(is.finite(alpha)) || any(alpha < 0) || sum(alpha) >= 1) {
            alpha <- NULL
            warning("invalid parameter 'alpha': trimming is not applied")
        }
    } else alpha <- NULL
    # observations with missings are excluded from simulation
    exclude <- getExclude(dataS)       
    if(length(exclude)) dataS <- dataS[-exclude,]
    # variables are coerced to factors
    select <- c(strata, basic)
    dataS <- checkFactor(dataS, select)
    dataP <- checkFactor(dataP, select)
    # sample data of variable to be simulated
    additionalS <- dataS[, additional]
    
    ## determine which models to fit and do further initializations
    haveBreaks <- !is.null(breaks)
    if(method == "multinom") {
        useMultinom <- TRUE
        useLogit <- FALSE
        useLm <- FALSE
        # define break points (if missing)
        if(haveBreaks) {
            checkBreaks(breaks)
            breaks <- if(zeros) union(breaks, 0) else unique(breaks)
            breaks <- sort(breaks)
        } else {
            if(is.null(upper) && gpd) upper <- Inf
            breaks <- getBreaks(additionalS, dataS[, w], 
                zeros, lower, upper, equidist, probs)
        }
    } else {
        useLm <- TRUE
        if(log) {
            if(is.null(const)) {
                ## use log-transformation
                # check for negative values
                neg <- which(additionalS < 0)
                haveNeg <- length(neg) > 0
                if(haveNeg) {
                    # define break points for negative values
                    if(haveBreaks) {
                        checkBreaks(breaks)
                        breaks <- c(unique(breaks[breaks < 0]), 0)
                    } else {
                        breaks <- getBreaks(additionalS[neg], dataS[neg, w], 
							zeros=TRUE, lower, upper)
                    }
                    if(zeros || length(breaks) > 2) {
                        useMultinom <- TRUE
                        breaks <- c(breaks, Inf)  # add Inf to breakpoints
                    } else useMultinom <- FALSE
                    useLogit <- !useMultinom
                } else {
                    useLogit <- zeros || any(additionalS == 0)
                    useMultinom <- FALSE
                }
            } else {
                # check constant
                if(!is.numeric(const) || length(const) == 0) {
                    stop("'const' must be numeric")
                } else const <- const[1]
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
    indStrata <- split(indP, dataP[, strata, drop=FALSE])
    # preparations for formulas and models
    predNames <- basic  # names of predictor variables
    fpred <- paste(predNames, collapse = " + ")  # for formula
	# check if population data contains factor levels that do not exist 
	# in the sample
	newLevels <- lapply(predNames, 
		function(nam) {
			levelsS <- levels(dataS[, nam])
			levelsP <- levels(dataP[, nam])
			levelsP[!(levelsP %in% levelsS)]
		})
	hasNewLevels <- sapply(newLevels, length) > 0
	excludeLevels <- any(hasNewLevels)
	
    ## preparations for multinomial or binomial logit model
    if(useMultinom || useLogit) {
        name <- getCatName(additional)
        fstring <- paste(name , "~" , fpred)  # formula for model as string
    }
    
    if(useMultinom) {
        ## some preparations
        dataS[, name] <- getCat(additionalS, breaks, zeros, right=TRUE)
        response <- dataS[, name]  # response variable
        # check threshold for GPD (if supplied)
        if(!useLm && gpd && !is.null(threshold) && length(threshold) != 1) {
            stop("'threshold' must be a single numeric value")
        }
        
        ## simulate categories
        # TODO: share code with 'simCategorical'
        valuesCat <- lapply(levels(dataS[, strata]), 
            function(s) {
				# sample data
				dataSample <- dataS[dataS[, strata] == s, , drop=FALSE]
				if(!nrow(dataSample)) return(character())
				# population data
				dataPop <- dataP[indStrata[[s]], predNames, drop=FALSE]
				# unique combinations in the stratum of the population need 
				# to be computed for prediction
				indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
				grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]
				# in sample, observations with NAs have been removed to fit the 
				# model, hence population can have additional levels
				# these need to be removed since those probabilities cannot 
				# be predicted from the model
				if(excludeLevels) {
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
				command <- paste("suppressWarnings(multinom(", fstring, 
					", weights=", w, ", data=dataSample, trace=FALSE", 
					", maxit=maxit, MaxNWts=MaxNWts))", sep="")
				mod <- eval(parse(text=command))  # fitted model
				# predict probabilities
				if(length(exclude) == 0) {
					probs <- predict(mod, newdata=grid, type="probs")
				} else {
					probs <- predict(mod, newdata=grid[-exclude, , drop=FALSE], 
						type="probs")
				}
                # set too small probabilities to exactly 0
                if(!is.null(eps)) {
                    probs[probs < eps] <- 0
                }
				# ensure it works for missing levels of response
				ind <- as.integer(which(table(dataSample[, name]) > 0))
				if(length(ind) > 2 && (nrow(grid)-length(exclude)) == 1) {
					probs <- t(probs)
				}
                # account for structural zeros
                if((!is.null(limit) || !is.null(censor)) && !is.null(dim(probs))) {
                    if(length(exclude) == 0) {
                        probs <- adjustProbs(probs, grid, names(indGrid), 
                            limit, censor)
                    } else {
                        probs <- adjustProbs(probs, grid[-exclude, , drop=FALSE], 
                            names(indGrid)[-exclude], limit, censor)
                    }
                }
				# local function for sampling from probabilities
				if(length(ind) == 1) {
					resample <- function(k, n, p) rep.int(1, n[k])
				} else if(length(ind) == 2) {
					resample <- function(k, n, p) spSample(n[k], c(1-p[k],p[k]))
				} else resample <- function(k, n, p) spSample(n[k], p[k,])
				# generate realizations for each combination
				sim <- as.list(rep.int(NA, length(indGrid)))
				if(length(exclude) == 0) {
					ncomb <- as.integer(sapply(indGrid, length))
					sim <- lapply(1:length(ncomb), resample, ncomb, probs)
				} else {
					ncomb <- as.integer(sapply(indGrid[-exclude], length))
					sim[-exclude] <- lapply(1:length(ncomb), 
						resample, ncomb, probs)
				}
				sim <- unsplit(sim, dataPop, drop=TRUE)
				# return realizations
				levels(response)[ind][sim]
            })
        valuesCat <- factor(unsplit(valuesCat, dataP[, strata, drop=FALSE]), 
            levels=levels(response))
        
        ## simulate (semi-)continuous values
        tcat <- table(valuesCat)
        ncat <- length(tcat)
        icat <- 1:ncat
        values <- as.list(rep.int(NA, ncat))
        # zeros
        if(zeros) {
            izero <- which(breaks == 0)
            values[izero] <- 0
            tcat <- tcat[-izero]
            ncat <- length(tcat)
            icat <- icat[-izero]
        }
        # values to be simulated with linear model or draws from Pareto 
        # distribution
        if(useLm) {
            # last breakpoint is Inf, the one before is 0
            nunif <- ncat - 1  # leave category of positive values out
        } else {
            nbreaks <- length(breaks)
            if(gpd) {
                if(is.null(threshold)) {
                    if(!haveBreaks && (!isTRUE(equidist) || !is.null(probs))) {
						ngpd <- nbreaks-2
					} else ngpd <- nbreaks-1
                } else if(any(tmp <- breaks >= threshold)) ngpd <- min(which(tmp))
                else ngpd <- nbreaks
            } else ngpd <- nbreaks
            if(gpd && ngpd <= ncat) {
                # adjust threshold and fit GPD
                threshold <- breaks[ngpd]  # adjust threshold
                estPar <- fitgpd(additionalS, threshold, est)  # fit GPD
                estPar <- estPar[["fitted.values"]]  # parameters of GPD
                # generalized pareto distribution
                igpd <- ngpd:ncat
                values[icat[igpd]] <- lapply(igpd, function(i) {
                        truncPareto(tcat[i], loc=threshold, 
                            scale=estPar["scale"], shape=estPar["shape"], 
                            breaks[i], breaks[i+1])
                    })
            }
            nunif <- ngpd - 1
        }
        # uniform distribution
        if(nunif > 0) {
            iunif <- 1:nunif
            values[icat[iunif]] <- lapply(iunif, function(i) {
                    runif(tcat[i], breaks[i], breaks[i+1])
                })
        }
        # turn list into vector of values
        values <- unsplit(values, valuesCat)
        
    } else if(useLogit) {
        ## some preparations
        if(log && is.null(const) && haveNeg) {
			indS <- additionalS > 0 
		} else indS <- additionalS != 0
        dataS[, name] <- as.integer(indS)
        formula <- as.formula(fstring)  # formula for model
        # auxiliary model for all strata (used in case of empty combinations)   
        useAux <- !is.null(tol)
        if(useAux) {
            if(length(tol) != 1 || tol <= 0) {
                stop("'tol' must be a single small positive value")
            }
			X <- model.matrix(formula, data=dataS)
			y <- dataS[, name]
            weights <- dataS[, w]
            mod <- logitreg(X, y, weights=weights)
            par <- mod$par
        }
        
        ## simulate binary vector
        valuesCat <- lapply(levels(dataS[, strata]), 
            function(s) {
                # sample data
                dataSample <- dataS[dataS[, strata] == s, , drop=FALSE]
                if(!nrow(dataSample)) return(numeric())
				# population data
				dataPop <- dataP[indStrata[[s]], predNames, drop=FALSE]
				# unique combinations in the stratum of the population need 
				# to be computed for prediction
				indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
				grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]
				# in sample, observations with NAs have been removed to fit the 
				# model, hence population can have additional levels
				# these need to be removed since those probabilities cannot 
				# be predicted from the model
				if(excludeLevels) {
					exclude <- mapply(function(pop, new) pop %in% new, 
						pop=grid[, hasNewLevels, drop=FALSE], 
						new=newLevels[hasNewLevels])
					if(is.null(dim(exclude))) {
						exclude <- which(any(exclude))
					} else exclude <- which(apply(exclude, 1, any))
					if(length(exclude) > 0) grid <- grid[exclude, , drop=FALSE]
					for(j in predNames[hasNewLevels]) {
						# drop new factor levels
						grid[, j] <- factor(as.character(grid[, j]), 
							levels=levels(dataS[, j]))
					}
				}
				# add 0 variable to combinations for use of 'model.matrix' 
				Xnew <- cbind(grid, 0)
				names(Xnew) <- c(predNames, name)
				Xnew <- model.matrix(formula, data=Xnew)
				# fit logit model
				X <- model.matrix(formula, data=dataSample)
				y <- dataSample[, name]
				weights <- dataSample[, w]
				mod <- logitreg(X, y, weights=weights)
				# add parameters from auxiliary model if necessary
				if(useAux) {
					indPar <- abs(mod$par) < tol
					mod$par[indPar] <- par[indPar]
				}
				# predict probabilities
				tmp <- exp(Xnew %*% mod$par)
				# avoid integer overflow
				p <- ifelse(is.infinite(tmp), 1, as.numeric(tmp / (1 + tmp)))
                # set too small probabilities to exactly 0
                if(!is.null(eps)) {
                    p[p < eps] <- 0
                }
				# generate realizations for each combination
				if(length(exclude) == 0) {
					ncomb <- as.integer(sapply(indGrid, length))
					sim <- lapply(1:length(ncomb), 
						function(k) spSample(ncomb[k], c(1-p[k], p[k])) - 1)
				} else {
					ncomb <- as.integer(sapply(indGrid[-exclude], length))
					sim <- as.list(rep.int(NA, length(indGrid)))
					sim[-exclude] <- lapply(1:length(ncomb), 
						function(k) spSample(ncomb[k], c(1-p[k], p[k])) - 1)
				}
				# return realizations
				unsplit(sim, dataPop, drop=TRUE)
            })
        valuesCat <- unsplit(valuesCat, dataP[, strata, drop=FALSE])
    }
    
    if(useLm) {
        ## some preparations
        if(useMultinom) {
            catLm <- names(tcat)[ncat]  # category for positive values
            dataS <- dataS[response == catLm, , drop=FALSE]
            indP <- valuesCat == catLm  
        } else if(useLogit) {
            dataS <- dataS[indS, , drop=FALSE]  # select only non-zeros
            indP <- valuesCat == 1  # indicates non-zeros in population
        }
        if(useMultinom || useLogit) {
            # adjust population data
            dataPop <- dataP[which(indP), , drop=FALSE]
            # list indStrata is adjusted so that it only contains
            # indices of persons in population with non-zero value
            indStrata <- split(1:nrow(dataPop), dataPop[, strata, drop=FALSE])
        } else dataPop <- dataP
        
        ## trim data (if specified)
        if(!is.null(alpha)) {
            additionalS <- dataS[, additional]
            p <- c(alpha[1], 1-alpha[2])
            bounds <- quantileWt(additionalS, dataS[, w], p)
            select <- additionalS > bounds[1] & additionalS < bounds[2]
            dataSample <- dataS[select, , drop=FALSE]
            
            # check if all relevant levels of predictor variables are still 
            # contained in sample after trimming
            # if not, trimming is not applied and a warning message is generated
            check <- unlist(sapply(predNames, function(i) {
                        table(dataS[, i]) > 0 & table(dataSample[, i]) == 0
                    }))
            if(any(check)) {
                dataSample <- dataS
                warning("trimming could not be applied")
            }
        } else dataSample <- dataS
        
        ## fit linear model
        # formula for linear model
        if(log) {
            fname <- paste("log(", additional, 
                if(!is.null(const)) " + const", ")", sep = "") 
        } else fname <- additional
        fstring <- paste(fname, " ~ ", fpred, sep = "")
        formula <- as.formula(fstring)
        # auxiliary model for all strata (used in case of empty combinations)
        weights <- dataSample[, w]
        mod <- lm(formula, weights=weights, data=dataSample)
        coef <- coef(mod)
        # simulate values
        valuesTmp <- lapply(levels(dataS[, strata]), 
            function(s, coef) {
                # sample data
                dataSample <- dataSample[dataSample[, strata] == s, , drop=FALSE]
                if(!nrow(dataSample)) return(numeric())
                # population data
				dataPop <- dataPop[indStrata[[s]], predNames, drop=FALSE]
				# unique combinations in the stratum of the population need 
				# to be computed for prediction
				indGrid <- split(1:nrow(dataPop), dataPop, drop=TRUE)
				grid <- dataPop[sapply(indGrid, function(i) i[1]), , drop=FALSE]
				# in sample, observations with NAs have been removed to fit the 
				# model, hence population can have additional levels
				# these need to be removed since those probabilities cannot 
				# be predicted from the model
				if(excludeLevels) {
					exclude <- mapply(function(pop, new) pop %in% new, 
						pop=grid[, hasNewLevels, drop=FALSE], 
						new=newLevels[hasNewLevels])
					if(is.null(dim(exclude))) {
						exclude <- which(any(exclude))
					} else exclude <- which(apply(exclude, 1, any))
					if(length(exclude) > 0) grid <- grid[exclude, , drop=FALSE]
					for(j in predNames[hasNewLevels]) {
						# drop new factor levels
						grid[, j] <- factor(as.character(grid[, j]), 
							levels=levels(dataS[, j]))
					}
				}
				# fit linear model
				#weights <- dataSample[, w]
				# command needs to be constructed as string
				command <- paste("lm(", fstring, 
					", weights=", w, ", data=dataSample)", sep="")
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
				if(length(exclude) == 0) {
					pred <- spPredict(mod, newdata)
				} else {
					pred <- as.list(rep.int(NA, length(indGrid)))
					pred[-exclude] <- spPredict(mod, newdata)
				}
				pred <- unsplit(pred, dataPop, drop=TRUE)
				# add error terms
                if(residuals) {
                    error <- sample(residuals(mod), 
                        size=nrow(dataPop), replace=TRUE)
                } else {
                    mu <- median(residuals(mod))
                    sigma <- mad(residuals(mod))
                    error <- rnorm(nrow(dataPop), mean=mu, sd=sigma)
                }
                # return realizations
                sim <- pred + error
                if(log) {
                    res <- exp(sim)  # transform back
                    if(!is.null(const)) res <- res - const  # subtract constant
                    res
                } else sim
            }, coef)
        valuesTmp <- unsplit(valuesTmp, dataPop[, strata, drop=FALSE])
        
        ## put simulated values together
        if(useMultinom) {
            values[which(indP == 1)] <- valuesTmp
        } else if(useLogit) {
            if(log && is.null(const) && haveNeg) {
                # only one category for non-positive values 
                # (two breakpoints, one of them is 0)
                values <- rep.int(NA, N)
                nonpos <- which(indP == 0)
                values[nonpos] <- runif(length(nonpos), breaks[1], breaks[2])
            } else values <- ifelse(is.na(indP), NA, 0) # only zeros
            values[which(indP == 1)] <- valuesTmp
        } else values <- valuesTmp
        
    }
    
    # attach new variable(s) to population data
    if(useMultinom && keep) dataP[, name] <- valuesCat
    dataP[, additional] <- values
    
    # return simulated data
    return(dataP)
}

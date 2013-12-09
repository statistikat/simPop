# -----------------------
# Authors: Andreas Alfons
#          K.U.Leuven
# -----------------------

# TODO: further generalization (no or more than one stratification variable)

simRelation <- function(dataS, dataP, hid = "hhid", 
		w = "weight", strata = "region", 
		basic = c("age", "sex", "hsize"), 
		relation = "relate", head = "head", 
		direct = NULL, 
		additional = c("nation", "ethnic", "religion"), 
		limit = NULL, censor = NULL, maxit = 500, 
        MaxNWts = 2000, eps = NULL, seed) {
	
	##### initializations
	if(!missing(seed)) set.seed(seed)  # set seed of random number generator
	if(length(strata) != 1) { 
		stop("currently 'strata' must specify exactly one column of 'data'")
	}
	varNames <- c(hid=hid, w=w, strata=strata, basic, 
		relation=relation, additional)
	
	# check data
	if(all(varNames %in% names(dataS))) {
		dataS <- dataS[, varNames]
	} else stop("undefined variables in the sample data")
	if(!all(c(strata, basic, relation) %in% names(dataP))) {
		stop("undefined variables in the population data")
	}
	
	# observations with missings are excluded from simulation
	exclude <- getExclude(dataS)
	if(length(exclude)) dataS <- dataS[-exclude,]
	
	# variables are coerced to factors
	dataS <- checkFactor(dataS, c(strata, basic, relation, additional))
	dataP <- checkFactor(dataP, c(strata, basic, relation))
    
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
	
	##### simulation of variables using a sequence of multinomial models
	
	# predictor variables
	
	predNames <- basic  # names of predictor variables
	
	for(i in additional) {
#		# for testing purposes
#		print(paste(Sys.time(), i, sep=": "))
		
		# components of multinomial model are specified
		levelsResponse <- levels(dataS[, i])
		formula <- paste(i, "~", paste(predNames, collapse=" + "))
		
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
		
		# generate values of new variable
		values <- lapply(levels(dataS[, strata]), 
			function(s) {
				# sample data
				dataSample <- dataS[dataS[, strata] == s, , drop=FALSE]
				if(!nrow(dataSample)) return(character())
				
				# first step: simulate category of household head
				# sample data
				indSampleHead <- which(dataSample[, relation] == head)
				dataSampleHead <- dataSample[indSampleHead, , drop=FALSE]
				# population data
				dataPop <- dataP[indStrata[[s]], predNames, drop=FALSE]
				# limit population data to household heads
				# this is not stored in a separate data.frame to save memory
				relPop <- dataP[indStrata[[s]], relation]
				indPopHead <- which(relPop == head)
				dataPop <- dataPop[indPopHead, , drop=FALSE]
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
				command <- paste("suppressWarnings(multinom(", formula, 
					", weights=", w, ", data=dataSampleHead, trace=FALSE", 
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
				# ensure code works for missing levels of response
				ind <- as.integer(which(table(dataSampleHead[, i]) > 0))
                if(length(ind) > 2 && (nrow(grid)-length(exclude)) == 1) {
                    probs <- t(probs)
                }
                # account for structural zeros
                if((!is.null(limit) || !is.null(censor)) && !is.null(dim(probs))) {
                    if(length(exclude) == 0) {
                        probs <- adjustProbs(probs, grid, names(indGrid), 
                            limit[[i]], censor[[i]])
                    } else {
                        probs <- adjustProbs(probs, grid[-exclude, , drop=FALSE], 
                            names(indGrid)[-exclude], limit[[i]], censor[[i]])
                    }
                }
                # local function for sampling from probabilities
                if(length(ind) == 1) {
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
					sim[-exclude] <- lapply(1:length(ncomb), 
						resample, ncomb, probs)
				}
				sim <- unsplit(sim, dataPop, drop=TRUE)
				sim <- factor(levelsResponse[ind][sim], levels=levelsResponse)
				
				# second step: assign category of household head to 
				# directly related household members
				# we actually assign the category of the household head to
				# all household members (because we need that variable 
				# anyway) and replace the values of non-directly related 
				# in the third step
				hidPop <- dataP[indStrata[[s]], hid]
				names(sim) <- hidPop[indPopHead]
				sim <- sim[as.character(hidPop)]
				
				# third step: simulate category of non-directly related 
				# household members using category of household head as 
				# additional predictor
				# FIXME: it might be a problem that we can have new levels 
				#        for the household head that do not exist in the 
				#        sample and are hence not represented in the model
				# get indices of non-directly related household members
				indSampleNonDirect <- which(!(dataSample[, relation] %in% c(head, direct)))
				indPopNonDirect <- which(!(relPop %in% c(head, direct)))
				if(length(indSampleNonDirect) > 0 && length(indPopNonDirect) > 0) {
					# create additional predictor in the sample
					add <- dataS[indSampleHead, i]	# these are still numeric for population data too
					hidSample <- dataSample[, hid]
					names(add) <- hidSample[indSampleHead]
					add <- add[as.character(hidSample)]
					iHead <- getHeadName(i)
					dataSample[, iHead] <- add
					# limit sample and population data to non-directly related household members
					dataSample <- dataSample[indSampleNonDirect, , drop=FALSE]
					dataPop <- dataP[indStrata[[s]], predNames, drop=FALSE]
					dataPop[, iHead] <- sim
					dataPop <- dataPop[indPopNonDirect, , drop=FALSE]
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
					formula <- paste(formula, iHead, sep = " + ")
					command <- paste("suppressWarnings(multinom(", formula, 
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
                    # ensure code works for missing levels of response
					ind <- as.integer(which(table(dataSample[, i]) > 0))
                    if(length(ind) > 2 && (nrow(grid)-length(exclude)) == 1) {
                        probs <- t(probs)
                    }
                    # account for structural zeros
                    if((!is.null(limit) || !is.null(censor)) && !is.null(dim(probs))) {
                        if(length(exclude) == 0) {
                            probs <- adjustProbs(probs, grid, names(indGrid), 
                                limit[[i]], censor[[i]])
                        } else {
                            probs <- adjustProbs(probs, grid[-exclude, , drop=FALSE], 
                                names(indGrid)[-exclude], limit[[i]], censor[[i]])
                        }
                    }
                    # local function for sampling from probabilities
                    if(length(ind) == 1) {
                        resample <- function(k, n, p) rep.int(1, n[k])
                    } else if(length(ind) == 2) {
                        resample <- function(k, n, p) spSample(n[k], c(1-p[k],p[k]))
                    } else resample <- function(k, n, p) spSample(n[k], p[k,])
					# generate realizations for each combination
					if(length(exclude) == 0) {
						ncomb <- as.integer(sapply(indGrid, length))
						simTmp <- lapply(1:length(ncomb), resample, ncomb, probs)
					} else {
						ncomb <- as.integer(sapply(indGrid[-exclude], length))
						simTmp <- as.list(rep.int(NA, length(indGrid)))
						simTmp[-exclude] <- lapply(1:length(ncomb), 
							resample, ncomb, probs)
					}
					simTmp <- unsplit(simTmp, dataPop, drop=TRUE)
					simTmp <- levelsResponse[ind][simTmp]
					sim[indPopNonDirect] <- simTmp
				}
				
				# return realizations
				sim
			})
		values <- unsplit(values, dataP[, strata, drop=FALSE])
		
		## add new categorical variable to data set
		dataP[, i] <- values
		predNames <- c(predNames, i)
	}
	
	# return simulated data
	dataP
}

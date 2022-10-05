simulateValues <- function(dataSample, dataPop, params) {
   # dS <<- copy(dataSample)
   # dP <<- copy(dataPop)
   # pp <<- copy(params)
   # stop()
  if(params$verbose){
    message("starting simulateValues...")
  }
  if (!nrow(dataSample)) {
    return(character())
  }
  predNames <- NULL
  meth <- params$method
  cur.var <- params$cur.var
  hasNewLevels <- params$hasNewLevels
  newLevels <- params$newLevels
  head <- params$head
  excludeLevels <- params$excludeLevels
  w <- params$w
  relation <- params$relation
  eps <- params$eps
  formula.cmd <- params$formula.cmd
  formula.cmd_nd <- params$formula.cmd_nd
  limit <- params$limit[[cur.var]]
  censor <- params$censor[[cur.var]]
  levelsResponse <- params$levelsResponse
  hid <- params$hid
  direct <- params$direct
  TF_head_samp <- dataSample[[relation]] %in% head
  TF_direct_samp <- dataSample[[relation]] %in% direct
  TF_head_pop <- dataPop[[relation]] %in% head
  TF_direct_pop <- dataPop[[relation]] %in% direct

  # first step: simulate category of household head
  # sample data
  indSampleHead <- which(TF_head_samp)
  dataSampleWork <- dataSample[indSampleHead, ]

  # limit population data to household heads
  # this is not stored in a separate data.frame to save memory
  indPopHead <- which(TF_head_pop)
  dataPopWork <- dataPop[indPopHead,]
  dataPopWork[[hid]] <- NULL
  # unique combinations in the stratum of the population need
  # to be computed for prediction
  indGrid <- split(1:nrow(dataPopWork), dataPopWork, drop=TRUE)
  grid <- dataPopWork[sapply(indGrid, function(i) i[1])]

  # in sample, observations with NAs have been removed to fit the
  # model, hence population can have additional levels
  # these need to be removed since those probabilities cannot
  # be predicted from the model
  if (excludeLevels) {
    exclude <- mapply(function(pop, new) pop %in% new, pop=grid[, hasNewLevels, drop=FALSE], new=newLevels[hasNewLevels])
    if (is.null(dim(exclude))) {
      exclude <- which(any(exclude))
    } else {
      exclude <- which(apply(exclude, 1, any))
    }
  } else {
    exclude <- integer()
  }
  if (length(exclude) == 0) {
    newdata <- copy(grid)
  } else {
    newdata <- copy(grid[-exclude])
  }
  # fit multinomial model
  # command needs to be constructed as string
  # this is actually a pretty ugly way of fitting the model
  if(length(unique(dataSampleWork[[cur.var]]))>1){
    mod <- eval(parse(text=formula.cmd))  # fitted model
    # predict probabilities
    
    ind <- match(colnames(newdata), colnames(dataSample))
    for (i in 1:length(ind)) {
      if (is.factor(unlist(newdata[, i, with = FALSE]))) {
        newdata[,colnames(newdata)[i] := factor(as.character(unlist(newdata[,colnames(newdata)[i],with=FALSE])),levels(dataSample[[ind[i]]]))]
      }
    }
    if (meth %in% "multinom") {
      probs <- predict(mod, newdata = newdata, type = "probs")
    } else if ( meth %in% c("ctree","cforest") ) {
      probs <- predict(mod, newdata = data.table(newdata), type = "prob")
      probs <- do.call("rbind", probs)
    } else if (meth %in% c("ranger")) {
      probs <- predict(mod, data = newdata, type="response")$predictions
    }
  }else{
    probs <- rep(1L,length(newdata))
  }
  

  

  
  # set too small probabilities to exactly 0
  if (!is.null(eps)) {
    probs[probs < eps] <- 0
  }
  # ensure code works for missing levels of response
  ind <- as.integer(which(table(dataSampleWork[[cur.var]]) > 0))
  if (length(ind) > 2 && (nrow(grid)-length(exclude)) == 1) {
    probs <- t(probs)
  }
  # account for structural zeros
  if ((!is.null(limit) || !is.null(censor)) && !is.null(dim(probs))) {
    if (length(exclude) == 0) {
      probs <- adjustProbs(probs, grid, names(indGrid), limit, censor)
    } else {
      probs <- adjustProbs(probs, grid[-exclude, , drop=FALSE], names(indGrid)[-exclude], limit, censor)
    }
  }
  # local function for sampling from probabilities
  if (length(ind) == 1 || is.null(dim(probs))) {
    resample <- function(k, n, p) rep.int(1, n[k])
  } else {
    resample <- function(k, n, p) spSample(n[k], p[k, ])
  }
  # generate realizations for each combination
  if (length(exclude) == 0) {
    ncomb <- as.integer(sapply(indGrid, length))
    sim <- lapply(1:length(ncomb), resample, ncomb, probs)
  } else {
    ncomb <- as.integer(sapply(indGrid[-exclude], length))
    sim <- as.list(rep.int(NA, length(indGrid)))
    sim[-exclude] <- lapply(1:length(ncomb), resample, ncomb, probs)
  }
  sim <- unsplit(sim, dataPopWork, drop = TRUE)
  sim <- factor(levelsResponse[ind][sim], levels = levelsResponse)
  simHead <- dataPop[indPopHead][,c(hid), with = FALSE]
  simHead[,c(cur.var):=sim]
  # second step: assign category of household head to
  # directly related household members
  # we actually assign the category of the household head to
  # all household members (because we need that variable
  # anyway) and replace the values of non-directly related
  # in the third step
  dataPopHID <- dataPop[[hid]]
  sim <- simHead[list(dataPopHID), , on=c(hid)][[cur.var]]
  rm(simHead)

  # third step: simulate category of non-directly related
  # household members using category of household head as
  # additional predictor
  # FIXME: it might be a problem that we can have new levels
  #        for the household head that do not exist in the
  #        sample and are hence not represented in the model
  # get indices of non-directly related household members
  indSampleNonDirect <- which(!(params$TF_head_samp|params$TF_direct_samp))
  indPopNonDirect <- which(!(params$TF_head_pop|params$TF_direct_pop))
  if ( length(indSampleNonDirect) > 0 && length(indPopNonDirect) > 0 ) {
    # create additional predictor in the sample
    add <- dataSample[[cur.var]][indSampleHead] # these are still numeric for population data too
    hidSample <- dataSample[[hid]]
    names(add) <- hidSample[indSampleHead]
    add <- add[as.character(hidSample)]
    iHead <- getHeadName(cur.var)
    dataSample[, iHead] <- add
    # limit sample and population data to non-directly related household members
    dataSampleWork <- dataSample[indSampleNonDirect,]
    dataPop[, iHead] <- sim
    dataPopWork <- dataPop[indPopNonDirect,predNames]
    # unique combinations in the stratum of the population need
    # to be computed for prediction
    indGrid <- split(1:nrow(dataPopWork), dataPopWork, drop=TRUE)
    grid <- dataPopWork[sapply(indGrid, function(i) i[1]), , drop=FALSE]
    grid <- as.data.frame(grid)
    # in sample, observations with NAs have been removed to fit the
    # model, hence population can have additional levels
    # these need to be removed since those probabilities cannot
    # be predicted from the model
    if ( excludeLevels ) {
      exclude <- mapply(function(pop, new) pop %in% new, pop=grid[, hasNewLevels, drop=FALSE], new=newLevels[hasNewLevels])
      if ( is.null(dim(exclude)) ) {
        exclude <- which(any(exclude))
      } else {
        exclude <- which(apply(exclude, 1, any))
      }
    } else {
      exclude <- integer()
    }
    if (length(exclude) == 0) {
      newdata <- copy(grid)
    } else {
      newdata <- copy(grid[-exclude])
    }
    # fit multinomial model
    # command needs to be constructed as string
    # this is actually a pretty ugly way of fitting the model
    if(length(unique(dataSampleWork[[cur.var]]))>1){
      mod <- eval(parse(text=formula.cmd))  # fitted model 
      # predict probabilities
      
      ind <- match(colnames(newdata), colnames(dataSample))
      for (i in 1:length(ind)) {
        if (is.factor(unlist(newdata[[i]]))) {
          newdata[,colnames(newdata)[i] := factor(as.character(unlist(newdata[,colnames(newdata)[i],with=FALSE])),levels(dataSample[[ind[i]]]))]
        }
      }
      
      if (meth %in% "multinom" ) {
        probs <- predict(mod, newdata = newdata, type = "probs")
      } else if ( meth %in% c("ctree","cforest") ) {
        probs <- predict(mod, newdata = data.table(newdata), type = "prob")
        probs <- do.call("rbind", probs)
        if (ncol(probs) == 2) {
          probs <- probs[, 2]
        }
      } else if (meth %in% c("ranger")) {
        probs <- predict(mod, data = newdata, type = "response")$predictions
        colnames(probs) <- mod$forest$levels
      }
    }else{
      probs <- rep(1L,length(newdata))
    }
    

    
    # set too small probabilities to exactly 0
    if (!is.null(eps)) {
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
    if ( length(exclude) == 0 ) {
      ncomb <- as.integer(sapply(indGrid, length))
      simTmp <- lapply(1:length(ncomb), resample, ncomb, probs)
    } else {
      ncomb <- as.integer(sapply(indGrid[-exclude], length))
      simTmp <- as.list(rep.int(NA, length(indGrid)))
      simTmp[-exclude] <- lapply(1:length(ncomb), resample, ncomb, probs)
    }
    simTmp <- unsplit(simTmp, dataPopWork, drop=TRUE)
    simTmp <- levelsResponse[ind][simTmp]
    sim[indPopNonDirect] <- simTmp
  }
  if(params$verbose){
    message("done.")
  }
  if(any(is.na(sim))){
    print(summary(sim))
    stop("NA in return object from simulatevalues")
  }

  # return realizations
  return(sim)
}

#' Simulate categorical variables of population data
#'
#' Simulate categorical variables of population data taking relationships
#' between household members into account. The household structure of the
#' population data needs to be simulated beforehand using
#' [simStructure()].
#'
#' The values of a new variable are simulated in three steps, where the second
#' step is optional. First, the values of the household heads are simulated
#' with multinomial log-linear models. Second, individuals directly related to
#' the corresponding household head (as specified by the argument
#' `direct`) inherit the value of the latter. Third, the values of the
#' remaining individuals are simulated with multinomial log-linear models in
#' which the value of the respective household head is used as an additional
#' predictor.
#'
#' The number of cpus are selected automatically in the following manner. The
#' number of cpus is equal the number of strata. However, if the number of cpus
#' is less than the number of strata, the number of cpus - 1 is used by
#' default.  This should be the best strategy, but the user can also overwrite
#' this decision.
#'
#' @name simRelation
#' @param simPopObj a `simPopObj` containing population and household
#' survey data as well as optionally margins in standardized format.
#' @param relation a character string specifying the columns of `dataS`
#' and `dataP`, respectively, that define the relationships between the
#' household members.
#' @param head a character string specifying the category of the variable given
#' by `relation` that identifies the household head.
#' @param direct a character string specifying categories of the variable given
#' by `relation`. Simulated individuals with those categories directly
#' inherit the values of the additional variables from the household head. The
#' default is `NULL` such that no individuals directly inherit value from
#' the household head.
#' @param additional a character vector specifying additional categorical
#' variables of `dataS` that should be simulated for the population data.
#' @param limit this can be used to account for structural zeros. If only one
#' additional variable is requested, a named list of lists should be supplied.
#' The names of the list components specify the predictor variables for which
#' to limit the possible outcomes of the response. For each predictor, a list
#' containing the possible outcomes of the response for each category of the
#' predictor can be supplied. The probabilities of other outcomes conditional
#' on combinations that contain the specified categories of the supplied
#' predictors are set to 0. If more than one additional variable is requested,
#' such a list of lists can be supplied for each variable as a component of yet
#' another list, with the component names specifying the respective variables.
#' @param censor this can be used to account for structural zeros. If only one
#' additional variable is requested, a named list of lists or
#' `data.frame`s should be supplied. The names of the list components
#' specify the categories that should be censored. For each of these
#' categories, a list or `data.frame` containing levels of the predictor
#' variables can be supplied. The probability of the specified categories is
#' set to 0 for the respective predictor levels. If more than one additional
#' variable is requested, such a list of lists or `data.frame`s can be
#' supplied for each variable as a component of yet another list, with the
#' component names specifying the respective variables.
#' @param maxit,MaxNWts control parameters to be passed to
#' [nnet::multinom()] and [nnet::nnet()]. See the help file
#' for [nnet::nnet()].
#' @param nr_cpus if specified, an integer number defining the number of cpus
#' that should be used for parallel processing.
#' @param eps a small positive numeric value, or `NULL` (the default). In
#' the former case, estimated probabilities smaller than this are assumed to
#' result from structural zeros and are set to exactly 0.
#' @param seed optional; an integer value to be used as the seed of the random
#' number generator, or an integer vector containing the state of the random
#' number generator to be restored.
#' @param regModel allows to specify the variables or model that is used when
#' simulating additional categorical variables. The following choices are

#' available if different from `NULL`.
#'
#' - "basic": only the basic household variables (generated with [simStructure()]
#' are used.
#' - "available": all available variables (that are common in the sample and
#' the synthetic population such as previously generated variables) excluding
#' id-variables, strata variables and household sizes are used for the
#' modeling. This parameter should be used with care because all factors are
#' automatically used as factors internally.
#' - formula-object: users may also specify a formula (class 'formula') that
#' will be used. Checks are performed that all required variables are available.
#' If parameter `regModel` is `NULL`, only basic household variables are used
#' in any case.
#' @param verbose set to `TRUE` if additional print output should be shown.
#' @param method a character string specifying the method to be used for
#' simulating the additional categorical variables. Accepted values are
#'
#' - "multinom": estimation of the conditional probabilities using
#' multinomial log-linear models and random draws from the resulting
#' distributions
#' - "ctree": for using Classification trees
#' - "cforest": for using random forest (implementation in package party)
#' - "ranger": for using random forest (implementation in package ranger)
#' @param by defining which variable to use as split up variable of the estimation. Defaults
#' to the strata variable.
#' @return An object of class [simPopObj-class] containing survey
#' data as well as the simulated population data including the categorical
#' variables specified by `additional`.
#' @note The basic household structure needs to be simulated beforehand with
#' the function [simStructure()].
#' @author Andreas Alfons and Bernhard Meindl
#' @export
#' @seealso [simStructure()], [simCategorical()],
#' [simContinuous()], [simComponents()]
#' @keywords datagen
#' @md
#' @examples
#' data(ghanaS) # load sample data
#' samp <- specifyInput(
#'   data = ghanaS,
#'   hhid = "hhid",
#'   strata = "region",
#'   weight = "weight"
#' )
#' ghanaP <- simStructure(
#'   data = samp,
#'   method = "direct",
#'   basicHHvars = c("age", "sex", "relate")
#' )
#' class(ghanaP)
#'
#' \donttest{
#' ## long computation time ...
#' ghanaP <- simRelation(
#'   simPopObj = ghanaP,
#'   relation = "relate",
#'   head = "head",
#'   additional = c("nation", "ethnic", "religion"), nr_cpus = 1
#' )
#' str(ghanaP)
#'}
simRelation <- function(simPopObj, relation = "relate", head = "head",
  direct = NULL, additional,
  limit = NULL, censor = NULL, maxit = 500, MaxNWts = 2000,
  eps = NULL, nr_cpus=NULL, seed = 1, regModel = NULL, verbose = FALSE,
  method=c("multinom", "ctree","cforest","ranger"),
  by = "strata") {


  V1 <- x <- newAdditionalVarible <- NULL
  method <- match.arg(method)
  
  # set seed of random number generator
  if (!missing(seed)) {
    set.seed(seed, "L'Ecuyer")  # set seed of random number generator
  }
  dataS <- sampleObj(simPopObj)
  dataP <- popObj(simPopObj)
  data_pop <- popData(simPopObj)
  data_sample <- sampleData(simPopObj)
  basic <- simPopObj@basicHHvars
  w <- dataS@weight
  hid <- dataS@hhid
  ## checking for household heads
  problematicHHsample <- data_sample[, any(get(relation) == head), by = c(hid)][V1 == FALSE][[hid]]
  problematicHHpop <- data_pop[, any(get(relation) == head), by = c(hid)][V1 == FALSE][[hid]]
  makeOneRandom <- function(x, head) {
    x[sample(1:length(x), 1)] <- head[1]
    return(x)
  }
  if (length(problematicHHsample) > 0) {
    warning("Sample:There are households without a head, a random head will be asigned.")
    setkeyv(data_sample, hid)
    data_sample[list(problematicHHsample), c(relation) := makeOneRandom(get(relation), head), by = c(hid)]
  }
  if (length(problematicHHpop) > 0) {
    warning("Population:There are households without a head, a random head will be asigned.")
    setkeyv(data_pop, hid)
    data_pop[list(problematicHHpop), c(relation) := makeOneRandom(get(relation), head), by = c(hid)]
  }

  # parameters for parallel computing
  if (by == "strata") {
    curStrata <- dataS@strata
  } else{
    curStrata <- by
  }
  if (!curStrata %in% colnames(data_sample)) {
    stop(curStrata, " is defined as by variable, but not in the sample data set.")
  }
  if (!curStrata %in% colnames(data_pop)) {
    stop(curStrata, " is defined as by variable, but not in the population data set.")
  }


  nr_strata <- length(levels(data_sample[[curStrata]]))
  pp <- parallelParameters(nr_cpus=nr_cpus, nr_strata=nr_strata)

  parallel <- pp$parallel
  nr_cores <- pp$nr_cores
  have_win <- pp$have_win
  rm(pp)
  if (verbose) {
    cat("Dimension of the population:\n")
    print(dim(data_pop))
    cat("Dimension of the sample:\n")
    print(dim(data_sample))
  }
  if (any(additional %in% colnames(data_pop))) {
    stop("variables already exist in the population!\n")
  }
  if ((length(regModel) == 1 | inherits(regModel, "formula") )&
      length(additional) > 1) {
    if (inherits(regModel, "formula")) {
      regModelL <- list()
      for (i in seq_along(additional)) {
        regModelL[[i]] <- regModel
      }
      regModel <- regModelL
    } else if (regModel %in% c("available", "basic")) {
      regModel <- rep(regModel, length(additional))
    }
  } else if (inherits(regModel, "formula")) {
    regModelL <- list()
    for (i in seq_along(additional)) {
      regModelL[[i]] <- regModel
    }
    regModel <- regModelL
  }
  if (is.null(regModel)) {
    regModel <- rep("basic", length(additional))
  }
  if (verbose) {
    message("Used model formulas:\n")
    print(regModel)
    message("------------------------------ \n")
  }

  if (regModel[[1]] == "basic") {
    varNames <- unique(c(
      hid = hid,
      w = w,
      curStrata,
      basic,
      relation = relation,
      additional,
      basic
    ))
  } else{
    varNames <- unique(c(
      hid = hid,
      w = w,
      curStrata,
      basic,
      relation = relation,
      additional,
      labels(terms(regModel[[1]]))
    ))
  }
  # check data
  if ( all(varNames %in% names(data_sample)) ) {
    data_sample <- data_sample[, varNames, with=F]
  } else {
    stop("undefined variables in the sample data\n")
  }
  if ( !all(c(curStrata, basic, relation) %in% names(data_pop)) ) {
    stop("undefined variables in the population data\n")
  }


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

  # list indStrata contains the indices of data_pop split by strata
  N <- nrow(data_pop)
  indStrata <- split(1:N, data_pop[[curStrata]])

  ##### simulation of variables using a sequence of multinomial models
  if( !missing(seed) ) {
    set.seed(seed,"L'Ecuyer")  # set seed of random number generator
  }

  # predictor variables
  counter <- 0
  for ( i in additional ) {
    counter <- counter+1
    if(verbose) cat(paste0("Simulating variable '",i,"'.\n"))
    if(length(regModel)>1){
      curRegModel <- regModel[counter]
    }else{
      curRegModel <- regModel
    }
    regInput <- regressionInput(simPopObj, additional=additional[counter], regModel=curRegModel)
    # names of predictor variables
    predNames <- setdiff(regInput[[1]]$predNames, c(dataS@hhsize, curStrata, relation))

    # observations with missings are excluded from simulation
    exclude <- getExclude.data.table(data_sample[,c(additional,predNames),with=FALSE])
    if ( length(exclude) > 0 ) {
      sampWork <- data_sample[-exclude,]
    } else {
      sampWork <- data_sample
    }

    # variables are coerced to factors
    sampWork <- checkFactor(sampWork, unique(c(curStrata, additional)))
    data_pop <- checkFactor(data_pop, unique(c(curStrata)))
    # components of multinomial model are specified
    levelsResponse <- levels(sampWork[[i]])
    # simulation of variables using a sequence of multinomial models
    if (method == "multinom") {
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd_nd <- paste(i, "~", paste(c(predNames,getHeadName(i)), collapse = " + "))
      formula.cmd <- paste0("suppressWarnings(multinom(", formula.cmd)
      formula.cmd_nd <- paste0("suppressWarnings(multinom(", formula.cmd_nd)
      if (!dataS@ispopulation) {
        formula.cmd <- paste0(formula.cmd,", weights=", dataS@weight)
        formula.cmd_nd <- paste0(formula.cmd_nd,", weights=", dataS@weight)
      }
      formula.cmd <- paste0(formula.cmd,", data=dataSampleWork, trace=FALSE",
                            ", maxit=",maxit, ", MaxNWts=", MaxNWts,"))")
      formula.cmd_nd <- paste0(formula.cmd_nd,", data=dataSampleWork, trace=FALSE",
                            ", maxit=",maxit, ", MaxNWts=", MaxNWts,"))")

      if (verbose) {
        cat("we are running the following multinom-model:\n")
        cat(strwrap(cat(gsub("))",")",gsub("suppressWarnings[(]","",formula.cmd)),"\n"), 76), sep = "\n")
      }
    } else if (method == "ctree") {
      # simulation via recursive partitioning and regression trees
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd_nd <- paste(i, "~", paste(c(predNames,getHeadName(i)), collapse = " + "))
      formula.cmd <- paste("suppressWarnings(ctree(", formula.cmd)
      formula.cmd_nd <- paste("suppressWarnings(ctree(", formula.cmd_nd)
      if (!dataS@ispopulation) {
        formula.cmd <- paste0(formula.cmd,", weights=as.integer(dataSampleWork$", dataS@weight,")")
        formula.cmd_nd <- paste0(formula.cmd_nd,", weights=as.integer(dataSampleWork$", dataS@weight,")")
      }
      formula.cmd <- paste0(formula.cmd, ", data=dataSampleWork))")
      formula.cmd_nd <- paste0(formula.cmd_nd, ", data=dataSampleWork))")
      if(verbose) {
        cat("we are running recursive partitioning:\n")
        cat(strwrap(cat(gsub("))",")",gsub("suppressWarnings[(]","",formula.cmd)),"\n"), 76), sep = "\n")
      }
    } else if (method == "cforest") {
      # simulation via random forest
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd_nd <- paste(i, "~", paste(c(predNames,getHeadName(i)), collapse = " + "))
      formula.cmd <- paste("suppressWarnings(cforest(", formula.cmd)
      formula.cmd_nd <- paste("suppressWarnings(cforest(", formula.cmd_nd)
      if (!dataS@ispopulation) {
        formula.cmd <- paste0(formula.cmd,", weights=as.integer(dataSampleWork$", dataS@weight,")")
        formula.cmd_nd <- paste0(formula.cmd_nd,", weights=as.integer(dataSampleWork$", dataS@weight,")")
      }
      formula.cmd <- paste0(formula.cmd,", data=dataSampleWork))")
      formula.cmd_nd <- paste0(formula.cmd_nd,", data=dataSampleWork))")
      if (verbose) {
        cat("we are running random forest classification (cforest):\n")
        cat(strwrap(cat(gsub("))",")",gsub("suppressWarnings[(]","",formula.cmd)),"\n"), 76), sep = "\n")
      }
    } else if (method == "ranger") {
      # simulation via random forest
      formula.cmd <- paste(i, "~", paste(predNames, collapse = " + "))
      formula.cmd_nd <- paste(i, "~", paste(c(predNames,getHeadName(i)), collapse = " + "))
      formula.cmd <- paste("suppressWarnings(ranger(", formula.cmd)
      formula.cmd_nd <- paste("suppressWarnings(ranger(", formula.cmd_nd)
      if (!dataS@ispopulation) {
        formula.cmd <- paste0(formula.cmd,", case.weights=dataSampleWork$", dataS@weight)
        formula.cmd_nd <- paste0(formula.cmd_nd,", case.weights=dataSampleWork$", dataS@weight)
      }
      formula.cmd <- paste0(formula.cmd, ", data=dataSampleWork,probability=TRUE))", sep="")
      formula.cmd_nd <- paste0(formula.cmd_nd, ", data=dataSampleWork,probability=TRUE))", sep="")
      if (verbose) {
        cat("we are running random forest (ranger):\n")
        cat(strwrap(cat(gsub("))",")",gsub("suppressWarnings[(]","",formula.cmd)),"\n"), 76), sep = "\n")
      }
    }
    newLevels <- lapply(predNames, function(nam) {
      levelsS <- levels(sampWork[[nam]])
      levelsP <- levels(data_pop[[nam]])
      levelsP[!(levelsP %in% levelsS)]
    })
    hasNewLevels <- sapply(newLevels, length) > 0
    excludeLevels <- any(hasNewLevels)

    params <- list()
    params$verbose <- verbose
    params$method <- method
    params$cur.var <- i
    params$hasNewLevels <- hasNewLevels
    params$newLevels <- newLevels
    params$strata <- dataS@strata
    params$curStrata <- curStrata
    params$head <- head
    params$excludeLevels <- excludeLevels
    params$w <- w
    params$relation <- relation
    params$formula.cmd <- formula.cmd
    params$formula.cmd_nd <- formula.cmd_nd
    params$eps <- eps
    params$limit <- limit
    params$censor <- censor
    params$levelsResponse <- levelsResponse
    params$hid <- hid
    params$direct <- direct
    if (parallel) {
      # windows
      if (have_win) {
        cl <- makePSOCKcluster(nr_cores)
        registerDoParallel(cl,cores=nr_cores)
        values <- foreach(x=levels(data_sample[[curStrata]]), .options.snow=list(preschedule=FALSE)) %dopar% {
          simulateValues(
            dataSample=data_sample[data_sample[[curStrata]]==x,],
            dataPop=data_pop[indStrata[[x]], c(predNames, simPopObj@pop@hhid, relation), with=FALSE], params
          )
        }
        stopCluster(cl)
      } else if (!have_win) {# linux/mac
        values <- mclapply(levels(data_sample[[curStrata]]), function(x) {
          simulateValues(
            dataSample=data_sample[data_sample[[curStrata]]==x,],
            dataPop=data_pop[indStrata[[x]], c(predNames, simPopObj@pop@hhid, relation), with=FALSE], params
          )
        }, mc.cores = min(nr_cores,length(levels(data_sample[[curStrata]]))))
      }
    } else {
      if (verbose) {
        message("Sample data:")
        summary(data_sample)
        message("Population data:")
        summary(data_pop)
      }

      values <- lapply(levels(data_sample[[curStrata]]), function(x) {
        simulateValues(
          dataSample=data_sample[list(x),,on=c(curStrata)],
          dataPop=data_pop[indStrata[[x]], c(predNames, simPopObj@pop@hhid, relation), with=FALSE], params
        )
      })
    }

    ## add new categorical variable to data set
    values <- factor(unsplit(values, data_pop[[curStrata]]), levels = levelsResponse)
    data_pop[[i]] <- values
  }
  # return simulated data
  simPopObj@pop@data <- data_pop
  invisible(simPopObj)
}

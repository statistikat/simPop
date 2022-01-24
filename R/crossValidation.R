# functionality similar to forcats::fct_unify
factor_unify <- function(fs){
  lev <- lapply(fs,levels)
  lev <- unique(do.call("c",lev))
  out <- lapply(fs,function(x,levNew){
    factor(as.character(x),levels=levNew,labels=levNew)
  },  levNew = lev)
  return(out)
}
categorical_metric <- function(x, y, weights) {
  unified_factors <- factor_unify(list(x = as.factor(x),
                                    y = as.factor(y)))
  x <- unified_factors$x
  y <- unified_factors$y
  Tx <- tableWt(x,
                weights = weights)
  Ty <- table(y)
  res <- mean((Tx - Ty)^2)
  return(res)
}


cross_validation <- function(synth_pop, fold = 1, grid,metric, sim, return_best = TRUE, verbose = FALSE) {
  
  fun_args <- as.list(formals(sim))
  fun_args$simPopObj <- synth_pop
  
  if(class(metric) != "function"){
    stop("metric is not a function, please provide a function")
  }
  
  first_level_params <- names(grid)[(names(grid) %in% names(fun_args))]
  second_level_params <- names(grid)[!(names(grid) %in% names(fun_args) | names(grid) == "fold")]
  
  weight <- synth_pop@sample@weight
  
  track <- cbind(fold_index = rep(1:fold, each = nrow(grid)),
                 do.call("rbind", replicate(fold, grid, simplify = FALSE)),
                 metric = NA)
  
  for (i in 1:nrow(track)) {
    
    for (first_level_param in first_level_params) {
      fun_args[first_level_param] <- track[i, first_level_param]
    }
    
    if(length(second_level_params) > 0){
      tmp_list <- as.list(track[i, second_level_params])
      names(tmp_list) <- second_level_params
      fun_args[["model_params"]] <- tmp_list
    }
    
    additional <- fun_args$additional
    
    ptm <- proc.time()
    synth_pop <- do.call(sim, fun_args)
    run_time <- proc.time() - ptm
    
    track[i, "run_time"] <- run_time[3]
    track[i, "metric"] <- metric(x = synth_pop@sample@data[,..additional][[1]],
                                 y = synth_pop@pop@data[,..additional][[1]], 
                                 weights = synth_pop@sample@data[, ..weight][[1]])
    
    if(verbose) {
      cat(paste0("    Metric at ", i," : ", round(track[i, "metric"], 4), "\n"))
    }
    
    
    dt_index <- names(synth_pop@pop@data)[grepl(names(synth_pop@pop@data), pattern = additional)]
    synth_pop@pop@data[, (dt_index) := NULL]
    
  }
  setDT(track)
  best <- track[,.(mean = mean(metric, na.rm = T), max = max(metric)), by=c(first_level_params, second_level_params)]
  setorder(best,mean)

  if(return_best){
    
    for (first_level_param in first_level_params) {
      fun_args[first_level_param] <- best[1, first_level_param]
    }
    
    if(length(second_level_params) > 0){
      tmp_list <- as.list(best[1, second_level_params])
      names(tmp_list) <- second_level_params
      fun_args[["model_params"]] <- tmp_list
    }
    
    synth_pop <- do.call(sim, fun_args)
    
    best_metric <- metric(x = synth_pop@sample@data[,..additional][[1]],
                          y = synth_pop@pop@data[,..additional][[1]], 
                          weights = synth_pop@sample@data[, ..weight][[1]])
    
    if(verbose) {
      cat(paste0("  Best metric : ", round(best_metric, 4), "\n"))
    }
    
    return(list(synth_pop = synth_pop, track = track))
  }else{
    return(track)
  }
}


#' Simulate variables of population data by cross validation
#'
#' Simulate variables of population data. The household structure
#' of the population data needs to be simulated beforehand.
#'
#' The number of cpus are selected automatically in the following manner. The
#' number of cpus is equal the number of strata. However, if the number of cpus
#' is less than the number of strata, the number of cpus - 1 is used by
#' default. This should be the best strategy, but the user can also overwrite
#' this decision.
#'
#' @name crossValidation
#' @param simPopObj a \code{simPopObj} containing population and household
#' survey data as well as optionally margins in standardized format.
#' @param additionals a character vector specifying additional categorical
#' variables available in the sample object of \code{simPopObj} that should be
#' simulated for the population data.
#' @param method a character string specifying the method to be used for
#' simulating the additional categorical variables. Accepted value at the moment
#' only 
#' \code{"xgboost"}  for using xgboost (implementation in package xgboost)
#' @param nr_cpus if specified, an integer number defining the number of cpus
#' that should be used for parallel processing.
#' @param regModel allows to specify the variables or model that is used when
#' simulating additional categorical variables. The following choices are
#' available if different from NULL.  \itemize{ \item'basic'only the basic
#' household variables (generated with \code{\link{simStructure}}) are used.
#' \item'available'all available variables (that are common in the sample and
#' the synthetic population such as previously generated varaibles) excluding
#' id-variables, strata variables and household sizes are used for the
#' modelling. This parameter should be used with care because all factors are
#' automatically used as factors internally.  \item formula-objectUsers may also
#' specify a specifiy formula (class 'formula') that will be used. Checks are
#' performed that all required variables are available.  } If method
#' 'distribution' is used, it is only possible to specify a vector of length
#' one containing one of the choices described above.  If parameter 'regModel'
#' is NULL, only basic household variables are used in any case.
#' @param verbose set to TRUE if additional print output should be shown.
#' @param by defining which variable to use as split up variable of the estimation. Defaults to the strata variable.
#' @param hyper_param_grid a grid which can contain model specific parameters which will be passed onto the function call for the respective model. 
#' @return An object of class \code{\linkS4class{simPopObj}} containing survey
#' data as well as the simulated population data including the categorical
#' variables specified by argument \code{additional}.
#' @note The basic household structure needs to be simulated beforehand with
#' the function \code{\link{simStructure}}.
#' @author Bernhard Meindl, Andreas Alfons, Stefan Kraft, Alexander Kowarik, Matthias Templ, Siro Fritzmann
#' @seealso \code{\link{simStructure}}, \code{\link{simRelation}},
#' \code{\link{simContinuous}}, \code{\link{simComponents}}, \code{\link{simCategorical}}
#' @keywords datagen
#' @examples
#' data(eusilcS) # load sample data
#' \donttest{
#' ## approx. 20 seconds computation time
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' ## in the following, nr_cpus are selected automatically
#' simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' grid <- expand.grid(nrounds = c(5, 10),
#'                     max_depth = 10,
#'                     eta = c(0.2, 0.3, 0.5),
#'                     eval_metric = "mlogloss",
#'                     stringsAsFactors = F)
#'
#' additionals <- c("pl030", "pb220a")
#' simPop <- crossValidation(simPop, additionals=c("pl030", "pb220a"), nr_cpus=1)
#' simPop
#' }
#' @export
crossValidation <- function(simPopObj, additionals, hyper_param_grid, fold = 3,
                            method = c("xgboost"), type = c("categorical"),
                            by = "strata", regModel = "available",
                            nr_cpus = 1, verbose = FALSE) {
  
  if(type == "categorical"){
    sim_pop_func <- simPop::simCategorical
    metric <- categorical_metric
    objective = "multi:softprob"
  } else {
    stop(paste0("Parameter type should be categorical, is ", type))
  }
  
  if(class(fold) == "numeric"){
    if(fold <= 0){
      stop(paste0("Parameter fold should be greater than 0, is ", fold))
    }
  } else{
    stop(paste0("Parameter fold should be class numeric, is ", class(fold)))
  }
  
  for (additional in additionals) {
    
    # TODO: only set if not exists, such that one could test multiple regModel
    # set additional, nr_cpus, regModel in grid
    hyper_param_grid$additional <- additional
    hyper_param_grid$nr_cpus <- nr_cpus
    hyper_param_grid$regModel <- regModel
    hyper_param_grid$by <- by
    hyper_param_grid$method <- method
    
    result_list <- cross_validation(simPopObj,
                                    grid = hyper_param_grid,
                                    fold = fold,
                                    metric = metric,
                                    sim = sim_pop_func,
                                    verbose = verbose,
                                    return_best = T)
    
    simPopObj <- result_list$synth_pop
    
    if(verbose){
      print(result_list$track)
    }
  }

  invisible(simPopObj)
}

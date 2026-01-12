#' Calibration of 0/1 weights by Simulated Annealing
#'
#' A Simulated Annealing Algorithm for calibration of synthetic population data
#' available in a \code{\linkS4class{simPopObj}}-object. The aims is to find,
#' given a population, a combination of different households which optimally
#' satisfy, in the sense of an acceptable error, a given table of specific
#' known marginals. The known marginals are also already available in slot
#' 'table' of the input object 'inp'.
#'
#' Calibrates data using simulated annealing. The algorithm searches for a
#' (near) optimal combination of different households, by swaping housholds at
#' random in each iteration of each temperature level. During the algorithm as
#' well as for the output the optimal (or so far best) combination will be
#' indicated by a logical vector containg only 0s (not inculded) and 1s
#' (included in optimal selection). The objective function for simulated
#' annealing is defined by the sum of absolute differences between target
#' marginals and synthetic marginals (=marginals of synthetic dataset). The sum
#' of target marginals can at most be as large as the sum of target marginals.
#' For every factor-level in \dQuote{split}, data must at least contain as many
#' entries of this kind as target marginals.
#'
#' Possible donors are automatically generated within the procedure.
#'
#' The number of cpus are selected automatically in the following manner. The
#' number of cpus is equal the number of strata. However, if the number of cpus
#' is less than the number of strata, the number of cpus - 1 is used by
#' default.  This should be the best strategy, but the user can also overwrite
#' this decision.
#'
#' @name calibPop
#' @docType data
#' @param data a \code{data.table} or \code{data.frame} holding the synthetic data or an object of class \code{\linkS4class{simPopObj}} where slot 'table' can optionally hold population constraints (see \code{\link{addKnownMargins}}).
#' @param hid column in \code{data} holding the household ID. If \code{NULL} a generic household ID will be created using \code{1:nrow(data)}. If \code{data} is an object of class \code{\linkS4class{simPopObj}} \code{hid} is extracted from the object.
#' @param pid column in \code{data} holding the personal ID. If \code{NULL} a generic personal ID will be created with each houeshold using \code{1:hsize}. If \code{data} is an object of class \code{\linkS4class{simPopObj}} \code{pid} is extracted from the object.
#' @param split given strata in which the problem will be split. Has to
#' correspond to a column population data (slot 'pop' of input argument 'inp')
#' . For example \code{split = (c("region")}, problem will be split for
#' different regions. Parallel computing is performed automatically, if
#' possible.
#' @param splitUpper optional column in the population for which decides the part
#' of the population from which to sample for each entry in \code{split}.
#' Has to correspond to a column population data (slot 'pop' of input argument 'inp').
#' For example \code{split = c("region"), splitUpper = c("Country")}
#' all units from the country are eligable for donor sample when problem is split
#' into regions. Is usefull if \code{simInitSpatial()} was used and the variable to split
#' the problem into results in very small groups (~couple of hundreds to thousands). 
#' @param temp starting temperatur for simulated annealing algorithm
#' @param epsP.factor a factor (between 0 and 1) specifying the acceptance
#' error for contingency table on individual level. For example epsP.factor = 0.05 results in an acceptance error for the
#' objective function of \code{0.05*sum(People)}.
#' @param epsH.factor a factor (between 0 and 1) specifying the acceptance
#' error for contingency table on household level. For example epsH.factor = 0.05 results in an acceptance error for the
#' objective function of \code{0.05*sum(Households)}.
#' @param epsMinN integer specifying the minimum number of units from which the synthetic populatin can deviate from cells in contingency tables.
#' This overwrites \code{epsP.factor} and \code{epsH.factor}. Is especially usefull if cells in \code{hhTables} and \code{persTables} are very small, e.g. <10.
#' @param maxiter maximum iterations during a temperature step.
#' @param temp.cooldown a factor (between 0 and 1) specifying the rate at which
#' temperature will be reduced in each step.
#' @param factor.cooldown a factor (between 0 and 1) specifying the rate at
#' which the number of permutations of housholds, in each iteration, will be
#' reduced in each step.
#' @param min.temp minimal temperature at which the algorithm will stop.
#' @param verbose boolean variable; if TRUE some additional verbose output is
#' provided, however only if \code{split} is NULL. Otherwise the computation is
#' performed in parallel and no useful output can be provided.
#' @param sizefactor the factor for inflating the population before applying 0/1 weights
#' @param choose.temp if TRUE \code{temp} will be rescaled according to \code{eps} and \code{choose.temp.factor}. \code{eps} is defined by the product between \code{epsP.factor} and \code{epsP.factor} with the sum over the target population margins supplied by \code{\link{addKnownMargins}} or \code{hhTables} and \code{persTables}. 
#' @param choose.temp.factor number between (0,1) for rescaling \code{temp} for simulated annealing. \code{temp} redefined by\code{max(temp,eps*choose.temp.factor)}.
#' Can be usefull if simulated annealing is split into subgroups with considerably different population sizes. Only used if \code{choose.temp=TRUE}.
#' @param scale.redraw Number between (0,1) scaling the number of households that need to be drawn and discarded in each iteration step.
#' The number of individuals currently selected through simulated annealing is substracted from the sum over the target population margins added to \code{inp} via \code{addKnownMargins}.
#' This difference is divided by the median household size resulting in an estimated number of housholds that the current synthetic population differs from the population margins (~\code{redraw_gap}).
#' The next iteration will then adjust the number of housholds to be drawn or discarded (\code{redraw}) according to \code{max(ceiling(redraw-redraw_gap*scale.redraw),1)} or \code{max(ceiling(redraw+redraw_gap*scale.redraw),1)} respectively.
#' This keeps the number of individuals in the synthetic population relatively stable regarding the population margins. Otherwise the synthetic population might be considerably larger or smaller then the population margins, through selection of many large or small households.
#' @param observe.times Number of times the new value of the objective function is saved. If \code{observe.times=0} values are not saved.
#' @param observe.break When objective value has been saved \code{observe.times}-times the coefficient of variation is calculated over saved values; if the coefficient of variation falls below \code{observe.break}
#' simmulated annealing terminates. This repeats for each new set of \code{observe.times} new values of the objecive function. Can help save run time if objective value does not improve much. Disable this termination by either setting \code{observe.times=0} or \code{observe.break=0}.
#' @param n.forceCooldown integer, if the solution does not move for \code{n.forceCooldown} iterations then a cooldown is automatically done.
#' @param hhTables Information on population margins for households. Can bei either a single \code{data.table} or \code{data.frame} or a list with multiple \code{data.tables}s or \code{data.frame}s.
#' Each table must have one column named \code{Freq} and all other columns holding variable(s) of the synthetic population. Each row in this table corresponds to a the frequency count a one of the variable combination in that table, see examples.
#' @param persTables Information on population margins for persons. Can bei either a single \code{data.table} or \code{data.frame} or a list with multiple \code{data.tables}s or \code{data.frame}s.
#' Each table must have one column named \code{Freq} and all other columns holding variable(s) of the synthetic population. Each row in this table corresponds to a the frequency count a one of the variable combination in that table, see examples.
#' @param redist.var single column in the population which can be redistributed in each `split`. Still experimental!
#' @param redist.var.factor numeric in the interval (0,1]. Used in combinationo with `redist.var`, still experimental!
#' @param ... arguments passed to \code{calibPop.data.table()}.
#' @return Returns an object of class \code{\linkS4class{simPopObj}} with an
#' updated population listed in slot 'pop'.
#' @author Bernhard Meindl, Johannes Gussenbauer and Matthias Templ
#' @references
#' M. Templ, B. Meindl, A. Kowarik, A. Alfons, O. Dupriez (2017) Simulation of Synthetic Populations for Survey Data Considering Auxiliary
#' Information. \emph{Journal of Statistical Survey}, \strong{79} (10), 1--38. \doi{10.18637/jss.v079.i10}
#' @keywords datasets
#'
#' @examples
#' data(eusilcS) # load sample data
#' data(eusilcP) # population data
#' \donttest{
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"), method="multinom", nr_cpus=1)
#'
#' # add margins
#' margins <- as.data.frame(
#'   xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$gender + eusilcP$citizenship))
#' colnames(margins) <- c("db040", "rb090", "pb220a", "Freq")
#' simPop <- addKnownMargins(simPop, margins)
#' }
#' 
#' epsP <- 0.1
#' # --------------------------------------------------
#' # apply simulated annealing using a simPop-Object
#' \donttest{
#' simPop_adj <- calibPop(simPop, split = "db040", temp = 1,
#'                        epsP.factor = epsP)
#' 
#' check_table <- merge(simPop_adj@pop@data[,.N,by=.(db040, rb090, pb220a)], 
#'                      margins, by=c("db040", "rb090", "pb220a"))
#' check_table[,sum(abs(N-Freq))/sum(N)<epsP]#' 
#' }
#' 
#' # --------------------------------------------------
#' # apply simulated annealing using a data.table
#' \donttest{
#' pop_data <- pop(simPop)
#' class(pop_data)
#' 
#' pop_adj <- calibPop(pop_data, hid = "db030", 
#'                    pid = "pid.simPop", split="db040", 
#'                     persTables = margins,  #  <- supply margins directly
#'                    temp=1, epsP.factor=epsP, choose.temp.factor = .5)
#' check_table <- merge(pop_adj[,.N,by=.(db040, rb090, pb220a)], 
#'                      margins, by=c("db040", "rb090", "pb220a"))
#' check_table[,sum(abs(N-Freq))/sum(N)<epsP]
#' 
#' }
#' 
#' # --------------------------------------------------
#' # apply simulated annealing using a data.frame
#' \donttest{
#' pop_data <- as.data.frame(pop(simPop))
#' class(pop_data)
#' 
#' pop_adj <- calibPop(pop_data, hid = "db030", pid = "pid.simPop", split="db040", 
#'                     persTables = margins,  #  <- supply margins directly
#'                     temp=1, epsP.factor=0.01, choose.temp.factor = .5)
#' margin_adj <- xtabs(rep(1, nrow(pop_adj)) ~ pop_adj$db040 + pop_adj$rb090 + pop_adj$pb220a)
#' margin_adj <- as.data.frame(margin_adj)
#' colnames(margin_adj) <- c("db040", "rb090", "pb220a", "N")
#' check_table <- merge(margin_adj, margins, by=c("db040", "rb090", "pb220a"))
#' sum(abs(check_table$N-check_table$Freq))/sum(check_table$N)<epsP
#' }
#' 
#' # --------------------------------------------------
#' # use multiple different margins
#' \donttest{
#' pop_data <- pop(simPop)
#' # person margins
#' persTables <- as.data.frame(
#'   xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$gender + eusilcP$citizenship))
#' colnames(persTables) <- c("db040", "rb090", "pb220a", "Freq")
#' 
#' # household margins
#' filter_hid <- !duplicated(eusilcP$hid)
#' eusilcP$hsize4 <- pmin(4,as.numeric(eusilcP$hsize))
#' hhTables <- as.data.frame(
#'   xtabs(rep(1, sum(filter_hid)) ~ eusilcP[filter_hid,]$region+eusilcP[filter_hid,]$hsize4))
#' colnames(hhTables) <- c("db040", "hsize4", "Freq")
#' pop_data[,hsize4 := pmin(4,as.numeric(hsize))]
#' 
#' epsP.factor <- 0.1
#' epsH.factor <- 0.1
#' 
#' pop_adj <- calibPop(pop_data, hid = "db030",
#'                     pid = "pid.simPop", split="db040",
#'                     temp=1, epsP.factor=0.1,
#'                     epsH.factor = 0.1,
#'                     persTables = persTables,
#'                     hhTables = hhTables)
#' 
#' check_pers <- merge(pop_adj[, .N, by = .(db040, rb090, pb220a)],
#'                     persTables, by = c("db040", "rb090", "pb220a"))
#' check_pers[, sum(abs(N-Freq))/sum(N) < epsP.factor]
#' 
#' check_house <- merge(pop_adj[!duplicated(db030), .N, by = .(db040, hsize4)],
#'                    hhTables, by = c("db040", "hsize4"))
#' check_house[, sum(abs(N-Freq))/sum(N) < epsH.factor]
#' }
#' 
#' @export
calibPop <- function(data, ...){
  UseMethod("calibPop")
}


#' @rdname calibPop
#' @export
calibPop.simPopObj <- function(data, ...){
  
  hid <- pid <- splitUpper <- zemp <- epsP.factor <- epsH.factor <- 
    epsMinN <- maxiter <- templ.cooldown <- factor.cooldown <- 
    min.temp <- sizefactor <- choose.temp <- choose.temp.factor <- 
    scale.redraw <- observe.times <- observe.break <- n.forceCooldown <-
    verbose <- hhTables <- persTables <- 
    temp <- temp.cooldown  <- redist.var <- redist.var.factor <- NULL
    
  
  
  # parse ellipsis
  ellipsis <- list(...)
  
  # get inputs from recordSwap.default
  rsArgs <- formals(calibPop.data.table)
  rsArgs$data <- rsArgs$`...` <- NULL
  
  # get parameters from ..., sdcObject and default values
  fun_params <- names(rsArgs)
  fun_values <- lapply(fun_params,function(z,ell,simPopObj,default){
    getVar(ell=ell,
           simPopObj=simPopObj,
           default=default,variable=z)
  },ell=ellipsis,simPopObj=data,default=rsArgs)
  names(fun_values) <- fun_params
  expr_values <- paste(fun_params,"<-",paste0("eval(fun_values[['",fun_params,"']])"))
  
  eval(parse(text = paste(expr_values)))

  # get data
  data0 <- popData(data)
  setDT(data0)
  
  # run record swaping default
  data0 <- calibPop.data.table(data = data0, hid = hid, pid = pid, 
                               split = split, splitUpper = splitUpper, 
                               temp = temp, epsP.factor = epsP.factor, epsH.factor = epsH.factor, 
                               epsMinN = epsMinN, maxiter = maxiter,
                               temp.cooldown = temp.cooldown, factor.cooldown = factor.cooldown, 
                               min.temp = min.temp, sizefactor = sizefactor, 
                               choose.temp = choose.temp, choose.temp.factor = choose.temp.factor,
                               scale.redraw = scale.redraw, observe.times = observe.times, observe.break = observe.break,
                               n.forceCooldown = n.forceCooldown, verbose = verbose,
                               hhTables = hhTables, persTables = persTables, redist.var = redist.var, 
                               redist.var.factor = redist.var.factor)
  
  
  data@pop@data <- data0
  invisible(data)
}

#' @rdname calibPop
#' @export
calibPop.data.frame <- function(data, ...){
  
  setDT(data)
  
  # run record swaping default
  data <- calibPop.data.table(data = data, ...)
  
  return(as.data.frame(data))
}

#' @rdname calibPop
#' @export
calibPop.data.table <- function(data, hid = attr(data, "hid"), pid = attr(data, "pid"),
                                split = NULL, splitUpper = NULL, temp = 1, epsP.factor = 0.05, epsH.factor = 0.05, 
                                epsMinN = 0, maxiter = 200, temp.cooldown = 0.9, factor.cooldown = 0.85, 
                                min.temp = 10^-3, sizefactor=2, choose.temp = TRUE, choose.temp.factor = 0.2,
                                scale.redraw = .5, observe.times = 50, observe.break = 0.05, n.forceCooldown = 100,
                                verbose = FALSE, hhTables = NULL, persTables = NULL, redist.var = NULL, redist.var.factor = 1, ...) {
  
  hid_simPop <- hhsize_simPop <- tmpsplit <- NULL
  
  if(verbose){
    t0 <- Sys.time()
  }
  
  # check hid
  if(is.null(hid)){
    data[,hid_simPop:=1:.N]
    hid <- "hid_simPop"
  }
  if(!hid %in% colnames(data)){
    stop("if data is of type 'data.table' or 'data.frame' then hid can either be 'NULL' or a column name in data!")      
  }
  data[,hhsize_simPop := .N,by=c(hid)]
  hhsize <- "hhsize_simPop"
  
  if ( length(split) > 1 ) {
    split <- split[1]
    warning("only first variable will be used to divide the population into strata")
  }
  if(is.null(splitUpper)){
    splitUpper <- split
  }
  if(length(splitUpper) > 1) {
    splitUpper <- splitUpper[1]
    warning("only first variable will be used to divide the population into strata")
  }
  
  weight_choose <- hid_help <- doub <- new.weights <- temporaryhid <- x <- NULL

  
  # ------------------------------------------------
  # check constraints

  if(is.null(persTables) & is.null(hhTables)){
    stop("Both persTables and hhTables are NULL.\nPopulation margins, either on individual or household level are needed to apply the calibration procedure!")
  }
  
  # check persTables
  if(!is.null(persTables)){
    if(verbose){
      cat("\nCheck Person-Tables\n")
    }
    persTables <- checkTables(persTables,data=data,split=split,verbose=verbose)
  }
  # check hhTables
  if(!is.null(hhTables)){
    if(verbose){
      cat("\nCheck Household-Tables\n")
    }
    hhTables <- checkTables(hhTables,data=data,split=split,verbose=verbose,namesTabs="hh")
  }
  totals <- c(persTables,hhTables)  
  totals <- totals[!sapply(totals,is.null)]
  # ------------------------------------------------
  
  # check some params
  if(!is.numeric(choose.temp.factor)){
    stop("choose.temp.factor must be numeric!")
  }
  choose.temp.factor <- choose.temp.factor[1]
  if(choose.temp.factor<=0|choose.temp.factor>=1){
    cat("choose.temp.factor must be in (0,1)\n Setting choose.temp.factor to default value of 0.2")
    choose.temp.factor <- .2
  }
  if(!is.numeric(scale.redraw)){
    stop("scale.redraw must be numeric!")
  }
  scale.redraw <- scale.redraw[1]
  if(scale.redraw<=0|scale.redraw>=1){
    cat("scale.redraw must be in (0,1)\n Setting scale.redraw to default value of 0.5")
    scale.redraw <- .5
  }
  
  if ( !is.null(split) ) {
    if ( !split %in% colnames(data) ) {
      stop("variable specified in argument 'split' must be a column in synthetic population (slot 'data' of argument 'inp')!\n")
    }
  } else {
    data$tmpsplit <- 1
    totals <- lapply(totals, function(z){
      z[,tmpsplit := 1]
    })
    split <- splitUpper <- "tmpsplit"
  }
  
  params <- list()
  params$temp <- as.numeric(temp)[1]
  params$epsP_factor = as.numeric(epsP.factor)[1]
  params$epsH_factor = as.numeric(epsH.factor)[1]
  params$maxiter = as.integer(maxiter)[1]
  params$temp_cooldown = as.numeric(temp.cooldown)[1]
  params$factor_cooldown = as.numeric(factor.cooldown)[1]
  params$min_temp = as.numeric(min.temp)[1]
  params$verbose <- ifelse(verbose, 1L, 0L)
  params$hhid <- hid
  params$pid <- pid
  params$hhsize <- hhsize
  params$epsMinN <- epsMinN
  params$redist.var <- redist.var
  params$redist.var.factor <- redist.var.factor
  if(!is.null(params[["redist.var"]])){
    params[["hhid_orig"]] <- paste0(params[["hhid"]],"_orig")
    params[["pid_orig"]] <- paste0(params[["pid"]],"_orig")
  }
  params$sizefactor <- sizefactor
  params$choose.temp <- choose.temp
  params$choose.temp.factor <- choose.temp.factor
  params$scale.redraw <- scale.redraw
  params$split <- split
  params$splitUpper <- splitUpper
  params$observe.times <- observe.times
  params$observe.break <- observe.break
  params$n.forceCooldown <- n.forceCooldown
  
  # make factor variables
  # and check if factor in totals coincide with factors in data
  data <- makeFactors(totals,data,split)
  
  # ## split the problem by "split"-factor
  # # parameters for parallel computing
  # strata <- unique(data[[split]])
  # nr_strata <- length(strata)
  # if(is.null(nr_cpus)){
  #   nr_cpus <- 1
  # }
  # nr_cpus <- min(nr_cpus, nr_strata, max(1,parallel::detectCores()-1))
  # setkeyv(data,split)
  # split_chunk <- paste0(split,"_simPop_chunk")
  # strata_on <- data.table(strata, rep(1:nr_cpus, length.out=nr_strata))
  # setnames(strata_on, colnames(strata_on), c(split, split_chunk))
  # data[strata_on, split_chunk := split_chunk, on=c(split), env = list(split_chunk = split_chunk)]
  # split.data <- split(data, f = data[[split_chunk]])

  
  split.number <- unique(data[[split]])
  final_weights <- lapply(1:length(split.number), function(x, data, totals, params) {
    split.level <- split.number[x]
    splitUpper.x <- data[list(split.level),,on=c(split)][[splitUpper]][1]
    data0 <- data[list(splitUpper.x),,on=c(splitUpper)]
    if(!is.null(params[["redist.var"]])){
      data0 <- multiply_data(data0,params,redist.var=params[["redist.var"]],split)
    }
    
    totals0 <- subsetList(totals,split=split,x=split.level)
    simAnnealingDT(
      data0 = data0,
      totals0 = totals0,
      split.level = split.level,
      params = params)
  }, data = data, totals = totals, params = params)

  # return dataset with new weights
  final_weights <- rbindlist(final_weights)
  final_weights <- final_weights[weight_choose>0]
  if(final_weights[,any(weight_choose>1)]){
    final_weights <- rbind(final_weights[weight_choose==1],final_weights[weight_choose>1,.SD[rep(1:.N,weight_choose)]])
  }
  final_weights[,doub:=1:.N,by=c(params[["pid"]],params[["hhid"]],split)]
  final_weights[,hid_help:=paste(get(params[["hhid"]]),get(split),doub,sep="_")]
  final_weights[,c(paste0(params[["hhid"]],"_new")):=.GRP,by=hid_help]
  final_weights[,c(paste0(params[["pid"]],"_new")):=paste(get(paste0(params[["hhid"]],"_new")),gsub("^[[:digit:]]*\\.","",get(params[["pid"]])),sep=".")]
  if(!is.null(params[["redist.var"]])){
    final_weights[,c(params[["hhid"]],params[["pid"]]):=NULL]
    setnames(final_weights,c(params[["hhid_orig"]],params[["pid_orig"]]),c(params[["hhid"]],params[["pid"]]))
  }
  
  final_weights[,c("weight_choose","doub","hid_help"):=NULL]
  
  data <- data[final_weights,,on=c(params[["hhid"]],params[["pid"]])]
  data[,c(params[["hhid"]],params[["pid"]],split,params[["redist.var"]]):=NULL]
  oldNames <- c(paste0(c(params[["hhid"]],params[["pid"]]),"_new"),paste0("i.",c(split,params[["redist.var"]])))
  newNames <- c(c(params[["hhid"]],params[["pid"]]),split,params[["redist.var"]])
  setnames(data,oldNames,newNames)
  
  if(verbose){
    dt <- difftime(Sys.time(),t0,units = "auto")
    cat("\nSimulated Annealing finished in: ",dt,attr(dt,"units"))
  }
  return(data)
}



#####
# additional helper functions
# subset list of data tables using on=.()
subsetList <- function(totals,split,x){
  
  totals <- lapply(totals,function(dt){
    dt.x <- dt[list(x),,on=c(split)]
    dt.x[[split]]<-NULL
    return(dt.x)
  })
  return(totals)
}  


# check list of contingency tables 
checkTables <- function(tabs,data,split=NULL,verbose=FALSE,namesTabs="pers"){
  
  if(is.null(tabs)){
    return(NULL)
  }
  
  if(!inherits(tabs, "list")){
    tabs <- list(tabs)
  }
  
  namesData <- copy(colnames(data))
  
  tabs <- lapply(tabs,function(z){
    
    if(!inherits(z, "data.frame")){
      stop("Each contingency Table must bei either a 'data.frame' or 'data.table'")  
    }
    setDT(z)
    
    if(!is.null(split) && !split%in%colnames(z)){
      stop(paste0("Variable ",split," not available in either persTables or hhTables!"))
    }
    
    if(!"Freq"%in%colnames(z)){
      stop(paste0("Totals, either number of households or number persons, in each contingency table must be specified by 'Freq'!"))
    }
    
    keepNames <- colnames(z)[colnames(z)%in%c(namesData,"Freq")]
    makeFactor <- keepNames[keepNames!="Freq"]
    z[,c(makeFactor):=lapply(.SD,factor),.SDcols=c(makeFactor)]
    
    z[, keepNames, with=FALSE]
    
  })
  
  names(tabs) <- paste0(namesTabs,1:length(tabs))
  
  if(verbose){
    cat("\nKeepig variables for Tables\n")
    sapply(tabs,function(z){
      print(colnames(z))
      cat("\n")
    })
  }
  
  if(verbose){
    cat("\nCheck if values in tables and data match\n")
  }
  lapply(tabs,function(z,data){
    check_values <- colnames(z)
    check_values <- check_values[check_values!="Freq"]
    tab_data <- data[,.N,by=c(check_values)]
    check_vars <- sapply(check_values,function(cz){
      isTRUE(all.equal(sort(as.character(unique(z[[cz]]))),sort(as.character(unique(tab_data[[cz]])))))
    })
    if(any(!check_vars)){
      stop("Values in columns ",paste(check_values,collapse=" and ")," do not agree between the synthetic population and the supplied population margins")
    }
  },data=data)
  
  
  
  return(tabs)
}

######################
# HELPERS
######
# make factor variables
# and check if factor in totals coincide with factors in dat
makeFactors <- function(totals,dat,split){
  ########
  # check and get all factors + levels in totals
  facLevels <- lapply(totals,function(z){
    
    facNames <- colnames(z)[which(sapply(z,is.factor))]
    output <- lapply(facNames,function(f){levels(z[[f]])})
    names(output) <- facNames
    return(output)
  })
  facLevels <- unlist(facLevels,recursive = FALSE)
  if(length(facLevels)==0){
    return(dat)
  }
  names(facLevels) <- gsub("^.*\\.","",names(facLevels))
  
  checkSplit <- facLevels[names(facLevels)==split]
  if(uniqueN(lengths(checkSplit))>1){
    stop(paste0("Splitvariable ",split," does not have the same number of values in each contingency table!"))
  }
  
  facLevels <- facLevels[!duplicated(facLevels)]
  
  ########
  # apply factor levels from totals to variables in dat
  # and make some checks
  
  for(i in 1:length(facLevels)){
    varName <- names(facLevels)[i]
    levels <- facLevels[[i]]
    dataLevels <- as.character(dat[,unique(get(varName))])
    
    leveldiff <- setdiff(dataLevels,levels)
    
    if(length(leveldiff)>0){
      stopMess <- paste0("Not all levels for variable ", varName," match in data and contingency table!\nLevels which only occure in data or contingency table:\n",paste(leveldiff,collapse=" "))
      stop(stopMess)
    }
    dat[,c(varName):=factor(get(varName),levels=levels)]
  }
  
  return(dat)
}

# multiply data if variable needs to be redistributed
multiply_data <- function(data0,params,redist.var,split){
  
  . <- hid <- NULL
  
  var_levels <- levels(data0[[redist.var]])
  fac_sample <- params[["redist.var.factor"]]
  redist.var_freq <- data0[,.(N_sample_extra=round(uniqueN(get(params[["hhid"]]))*fac_sample)),by=c(redist.var)]
  data0_extra <- list()
  for(u in 1:nrow(redist.var_freq)){
    uc <- redist.var_freq[u,][[redist.var]]
    N_samp <- redist.var_freq[u,][["N_sample_extra"]]
    
    samp_hid <- data0[!.(uc),unique(hid),on=c(redist.var)]
    samp_hid <- sample(samp_hid,min(length(samp_hid),N_samp),replace=FALSE)
    data0_u <- data0[.(samp_hid),on=.(hid)]
    set(data0_u,j=redist.var,value=uc)
    data0_extra <- c(data0_extra,list(data0_u))
  }
  data0_extra <- rbindlist(data0_extra)
  data0_extra[,c(split):=NA]
  data0 <- rbindlist(list(data0,data0_extra),use.names=TRUE)
  data0[,c(redist.var):=factor(get(redist.var),levels=var_levels)]
  
  var_values <- as.character(unique(data0[[redist.var]]))
  
  hid_orig <- params[["hhid_orig"]]
  pid_orig <- params[["pid_orig"]]
  setnames(data0,c(params[["hhid"]],params[["pid"]]),c(hid_orig,pid_orig))
  data0[,c(params[["hhid"]]):=.GRP,by=c(hid_orig,redist.var)]
  data0[,c(params[["pid"]]):=1:.N,by=c(params[["hhid"]])]
  
  return(data0)
}


# helpfunctino to get paramete values from ..., simPopObj and default values
getVar <- function(ell,simPopObj,default,variable){
  
  in_ell <- ell[[variable]]
  if(variable=="hid"){
    in_simPopObj <- popObj(simPopObj)@hhid
  }else if(variable=="pid"){
    in_simPopObj <- popObj(simPopObj)@pid
  }else if(variable=="persTables"){
    in_simPopObj <- tableObj(simPopObj)
    if("N" %in% colnames(in_simPopObj)){
      setnames(in_simPopObj,"N","Freq")
    }
  }else{
    in_simPopObj <- NULL
  }
  
  null_ell <- is.null(in_ell)
  null_simPopObj <- is.null(in_simPopObj)
  
  if(null_ell & null_simPopObj){
    # if(is.symbol(default[[variable]])){
    #   if(variable=="risk_variables"){
    #     stop("argument `",variable,"` is missing, with no default\n Alternatively one can specifcy `",variable,"` through the parameter `keyVars` in `createSdcObj()`")
    #   }else if(variable=="hid"){
    #     stop("argument `",variable,"` is missing, with no default\n Alternatively one can specifcy `",variable,"` through the parameter `hhId` in `createSdcObj()`")
    #   }else{
    #     stop("argument `",variable,"` is missing, with no default\n Alternatively one can specifcy `",variable,"` through the parameter `options` in `createSdcObj()`")
    #   }
    # }else{
    #   # set default value for variable
    #   take_value <- default[[variable]]
    # }
    take_value <- default[[variable]]
  }
  
  if(!null_ell & !null_simPopObj){
    warning("argument `",variable,"` defined in function call and in `data`: taking value from function call")
    take_value <- in_ell
  }
  
  if(!null_ell){
    take_value <- in_ell
  }else if(!null_simPopObj){
    take_value <- in_simPopObj
  }
  return(take_value)
}


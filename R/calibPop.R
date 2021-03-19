# required for parallelisation
calcFinalWeights <- function(data0, totals0, params) {
  # TODO: speed-improvements
  # each row of output contains indices of a specific margin
  indices_by_constraint <- function(data0, totals0, parameter) {
    # recode to character and set NA to "NA"
    myfun <- function(x) {
      x <- as.character(x)
      ii <- which(is.na(x))
      if ( length(ii) > 0 ) {
        x[ii] <- "NA"
      }
      x
    }

    out <- matrix(NA, nrow=nrow(totals0), ncol=nrow(data0))
    totals0 <- as.data.frame(totals0)
    totals0[,parameter] <- apply(totals0[,parameter,drop=FALSE], 2, myfun)
    for ( x in parameter ) {
      data0[[x]] <- myfun(data0[[x]])
    }
    for (i in 1:nrow(totals0) ) {
      ex <- paste("out[i,] <- as.integer(", sep="")
      for ( k in seq_along(parameter) ) {
        if ( k > 1 ) {
          ex <- paste(ex, "&", sep="")
        }
        ex <- paste(ex, 'data0[,"',parameter[k],'", with=FALSE] =="',totals0[i, parameter[k]],'"', sep="")
      }
      ex <- paste(ex, ")", sep="")
      eval(parse(text=ex))
    }
    out
  }

  inp <- indices_by_constraint(data0, totals0, params$parameter)
  current_totals <- as.numeric(totals0$N)
  weights <- as.integer(data0[[params$weight]])

  hh_info <- list()
  hh_info$hh_ids <- as.integer(data0[[params$hhid]])
  hh_info$hh_head <- rep(0L, nrow(data0))
  index <- which(sapply(data0[[params$pid]], function(x) {
    unlist(strsplit(x, "[.]"))[2] }) == "1")
  hh_info$hh_head[index] <- 1L
  hh_info$hh_size <- as.integer(data0[[params$hhsize]])
  hh_info$median_hhsize <- median(hh_info$hh_size[hh_info$hh_head==1], na.rm=TRUE)

  w <- calibPop_work(
    inp = inp,
    totals = current_totals,
    weights = weights,
    hh_info = hh_info,
    params = params
  )
  invisible(w)
}

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
#' @param inp an object of class \code{\linkS4class{simPopObj}} with slot
#' 'table' being non-null! (see \code{\link{addKnownMargins}}).
#' @param split given strata in which the problem will be split. Has to
#' correspond to a column population data (slot 'pop' of input argument 'inp')
#' . For example \code{split = c("region")}, problem will be split for
#' different regions. Parallel computing is performed automatically, if
#' possible.
#' @param temp starting temperatur for simulated annealing algorithm
#' @param eps.factor a factor (between 0 and 1) specifying the acceptance
#' error. For example eps.factor = 0.05 results in an acceptance error for the
#' objective function of \code{0.05*sum(totals)}.
#' @param maxiter maximum iterations during a temperature step.
#' @param temp.cooldown a factor (between 0 and 1) specifying the rate at which
#' temperature will be reduced in each step.
#' @param factor.cooldown a factor (between 0 and 1) specifying the rate at
#' which the number of permutations of housholds, in each iteration, will be
#' reduced in each step.
#' @param min.temp minimal temperature at which the algorithm will stop.
#' @param nr_cpus if specified, an integer number defining the number of cpus
#' that should be used for parallel processing.
#' @param verbose boolean variable; if TRUE some additional verbose output is
#' provided, however only if \code{split} is NULL. Otherwise the computation is
#' performed in parallel and no useful output can be provided.
#' @param sizefactor the factor for inflating the population before applying 0/1 weights
#' @param memory if TRUE simulated annealing is applied in less memory intensive way. Is especially usefull if factor or population is large. For this option simulated annealing is not entirely implemented in C++, therefore it might be slower than option \code{memory=FALSE}.
#' @param choose.temp if TRUE \code{temp} will be rescaled according to \code{eps} and \code{choose.temp.factor}. \code{eps} is defined by the product between \code{eps_factore} and the sum over the target population margins, see \code{\link{addKnownMargins}}. Only used if \code{memory=TRUE}.
#' @param choose.temp.factor number between (0,1) for rescaling \code{temp} for simulated annealing. \code{temp} redefined by\code{max(temp,eps*choose.temp.factor)}.
#' Can be usefull if simulated annealing is split into subgroups with considerably different population sizes. Only used if \code{choose.temp=TRUE} and \code{memory=TRUE}.
#' @param scale.redraw Only used if \code{memory=TRUE}. Number between (0,1) scaling the number of households that need to be drawn and discarded in each iteration step.
#' The number of individuals currently selected through simulated annealing is substracted from the sum over the target population margins added to \code{inp} via \code{addKnownMargins}.
#' This difference is divided by the median household size resulting in an estimated number of housholds that the current synthetic population differs from the population margins (~\code{redraw_gap}).
#' The next iteration will then adjust the number of housholds to be drawn or discarded (\code{redraw}) according to \code{max(ceiling(redraw-redraw_gap*scale.redraw),1)} or \code{max(ceiling(redraw+redraw_gap*scale.redraw),1)} respectively.
#' This keeps the number of individuals in the synthetic population relatively stable regarding the population margins. Otherwise the synthetic population might be considerably larger or smaller then the population margins, through selection of many large or small households.
#' @param observe.times Only used if \code{memory=TRUE}. Number of times the new value of the objective function is saved. If \code{observe.times=0} values are not saved.
#' @param observe.break Only used if \code{memory=TRUE}. When objective value has been saved \code{observe.times}-times the coefficient of variation is calculated over saved values; if the coefficient of variation falls below \code{observe.break}
#' simmulated annealing terminates. This repeats for each new set of \code{observe.times} new values of the objecive function. Can help save run time if objective value does not improve much. Disable this termination by either setting \code{observe.times=0} or \code{observe.break=0}.
#' @return Returns an object of class \code{\linkS4class{simPopObj}} with an
#' updated population listed in slot 'pop'.
#' @author Bernhard Meindl, Johannes Gussenbauer and Matthias Templ
#' @references
#' M. Templ, B. Meindl, A. Kowarik, A. Alfons, O. Dupriez (2017) Simulation of Synthetic Populations for Survey Data Considering Auxiliary
#' Information. \emph{Journal of Statistical Survey}, \strong{79} (10), 1--38. \doi{10.18637/jss.v079.i10}
#' @keywords datasets
#' @export
#' @examples
#' data(eusilcS) # load sample data
#' data(eusilcP) # population data
#' \donttest{
#' ## approx. 20 seconds computation time
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"), method="multinom", nr_cpus=1)
#'
#' # add margins
#' margins <- as.data.frame(
#'   xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$gender + eusilcP$citizenship))
#' colnames(margins) <- c("db040", "rb090", "pb220a", "freq")
#' simPop <- addKnownMargins(simPop, margins)
#' simPop_adj2 <- calibPop(simPop, split="db040", 
#'   temp=1, eps.factor=0.1,
#'   memory=TRUE, nr_cpus = 1)
#' }
#' # apply simulated annealing
#' \donttest{
#' ## long computation time
#' simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=FALSE)
#' }
calibPop <- function(inp, split, temp = 1, eps.factor = 0.05, maxiter=200,
  temp.cooldown = 0.9, factor.cooldown = 0.85, min.temp = 10^-3,
  nr_cpus=NULL, sizefactor=2, memory=TRUE,
  choose.temp=TRUE,choose.temp.factor=0.2,scale.redraw=.5,observe.times=50,observe.break=0.05,
  verbose=FALSE) {
  if(verbose){
    t0 <- Sys.time()
  }
  if ( class(inp) != "simPopObj" ) {
    stop("argument 'inp' must be of class 'simPopObj'!\n")
  }

  if ( length(split) > 1 ) {
    split <- split[1]
    warning("only first variable will be used to divide the population into strata")
  }

  hid_help <- doub <- new.weights <- temporaryhid <- x <- NULL
  data <- popData(inp)
  hid <- popObj(inp)@hhid
  pid <- popObj(inp)@pid
  hhsize <- popObj(inp)@hhsize
  totals <- tableObj(inp)
  parameter <- colnames(totals)[-ncol(totals)]

  if ( !is.null(split) ) {
    verbose <- FALSE
    if ( !split %in% colnames(data) ) {
      stop("variable specified in argument 'split' must be a column in synthetic population (slot 'data' of argument 'inp')!\n")
    }
    if ( !split %in% parameter ) {
      stop("variable specified in argument 'split' must be a column in slot 'table' of argument 'inp'!\n")
    }
  } else {
    data$tmpsplit <- 1
    totals$tmpsplit <- 1
    split <- "tmpsplit"
  }

  params <- list()
  params$temp <- as.numeric(temp)[1]
  params$eps_factor = as.numeric(eps.factor)[1]
  params$maxiter = as.integer(maxiter)[1]
  params$temp_cooldown = as.numeric(temp.cooldown)[1]
  params$factor_cooldown = as.numeric(factor.cooldown)[1]
  params$min_temp = as.numeric(min.temp)[1]
  params$verbose <- ifelse(verbose, 1L, 0L)
  params$parameter <- parameter
  params$weight <- popObj(inp)@weight
  params$hhid <- hid
  params$pid <- pid
  params$hhsize <- hhsize

  if(!memory){
    # generate donors
    data2 <- sampHH(data, sizefactor=sizefactor, hid=hid, strata=split, hsize=hhsize)
    data2[,(params$weight):=list(0)]
    setnames(data2,hid,"temporaryhid")
    if(is.numeric(data2[[hid]])){
      data2[,hid:=list(data[,sapply(.SD,max),.SDcols=hid]+temporaryhid)]
    }else{
      data2[,hid:=list(paste0("9999_",temporaryhid))]
    }
    setnames(data2,"temporaryhid",hid)
    #  data2 <- data2[, which(!grepl(hid, colnames(data2))), with=FALSE]
    cn <- colnames(data2)
    #  cn[length(cn)] <- hid
    #  setnames(data2, cn)
    data2 <- data2[,colnames(data), with=FALSE]
    #data2 <- data2[,match(colnames(data), cn), with=FALSE]
    data <- rbind(data, data2)
  }else{
    # check params for memory=TRUE
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
  }


  # parameters for parallel computing
  nr_strata <- length(unique(data[[split]]))
  pp <- parallelParameters(nr_cpus=nr_cpus, nr_strata=nr_strata)
  parallel <- pp$parallel
  nr_cores <- pp$nr_cores
  have_win <- pp$have_win; rm(pp)

  ii <- match(parameter, colnames(data))
  data[,ii] <- data[,lapply(.SD,as.factor),.SDcols=parameter]
  rm(ii)

  ## split problem by "split"-factor
  setkeyv(data, split)
  setkeyv(totals, split)
  split.number <- unique(data[,split,with=FALSE])

  if ( parallel ) {
    # windows
    if ( have_win ) {
      cl <- makePSOCKcluster(nr_cores)
      registerDoParallel(cl,cores=nr_cores)
      if(memory){
        final_weights <- foreach(x=1:nrow(split.number), .options.snow=list(preschedule=TRUE)) %dopar% {
          simAnnealingDT(
            data0=data[split.number[x]],
            totals0=totals[which(totals[,split,with=FALSE]==as.character(split.number[x][[split]])),],
            params=params,sizefactor=sizefactor,choose.temp=choose.temp,
            choose.temp.factor=choose.temp.factor,scale.redraw=scale.redraw,
            split.level=paste0(unlist(split.number[x])),observe.times=observe.times,observe.break=observe.break)
        }
      }else{
        final_weights <- foreach(x=1:nrow(split.number), .options.snow=list(preschedule=TRUE)) %dopar% {
          calcFinalWeights(
            data0=data[split.number[x]],
            totals0=totals[which(totals[,split,with=FALSE]==as.character(split.number[x][[split]])),],
            params=params
          )
        }
      }
      stopCluster(cl)
    }else if ( !have_win ) {# linux/mac
      if(memory){
        final_weights <- mclapply(1:nrow(split.number), function(x) {
          simAnnealingDT(
            data0=data[split.number[x]],
            totals0=totals[which(totals[,split,with=FALSE]==as.character(split.number[x][[split]])),],
            params=params,sizefactor=sizefactor,choose.temp=choose.temp,
            choose.temp.factor=choose.temp.factor,scale.redraw=scale.redraw,
            split.level=paste0(unlist(split.number[x])),observe.times=observe.times,observe.break=observe.break)
        },mc.cores=nr_cores)
      }else{
        final_weights <- mclapply(1:nrow(split.number), function(x) {
          calcFinalWeights(
            data0=data[split.number[x]],
            totals0=totals[which(totals[,split,with=FALSE]==as.character(split.number[x][[split]])),],
            params=params)
        },mc.cores=nr_cores)
      }
    }
  } else {
    if(memory){
      # use memory efficient but slower method
      final_weights <- lapply(1:nrow(split.number), function(x) {
        out <- simAnnealingDT(
          data0=data[split.number[x]],
          totals0=totals[which(totals[,split,with=FALSE]==as.character(split.number[x][[split]])),],
          params=params,sizefactor=sizefactor,choose.temp=choose.temp,
          choose.temp.factor=choose.temp.factor,scale.redraw=scale.redraw,
          split.level=paste0(unlist(split.number[x])),observe.times=observe.times,observe.break=observe.break)
        return(out)
      })
    }else{
      # use c++ implementation - can be quite memory intensive
      final_weights <- lapply(1:nrow(split.number), function(x) {
        calcFinalWeights(
          data0=data[split.number[x]],
          totals0=totals[which(totals[,split,with=FALSE]==as.character(split.number[x][[split]])),],
          params=params)
      })
    }
  }

  # return dataset with new weights
  data[,new.weights:=as.integer(unlist(final_weights))]
  if(memory){
    data <- data[new.weights>0,]
    data <- rbind(data[new.weights==1],data[new.weights>1,.SD[rep(1:.N,new.weights)]])
    data[,doub:=1:.N,by=c(params[["pid"]],params[["hhid"]])]
    data[,hid_help:=paste(get(params[["hhid"]]),doub,sep="_")]
    data[,c(params[["hhid"]]):=.GRP,by=hid_help]
    data[,c(params[["pid"]]):=paste(get(params[["hhid"]]),gsub("^[[:digit:]]*.","",get(params[["pid"]])),sep=".")]
    data[,c("new.weights","doub","hid_help"):=NULL]
  }else{
    data <- data[data$new.weights==1,]
    data[[inp@pop@weight]] <- data$new.weights
    data$new.weights <- NULL
  }

  inp@pop@data <- data
  if(verbose){
    Sys.time()-t0
  }
  invisible(inp)
}
NULL


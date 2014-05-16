# required for parallelisation
calcFinalWeights <- function(data0, totals0, params) {
  # TODO: speed-improvements
  # each row of output contains indices of a specific margin
  indices_by_constraint <- function(data0, totals0, parameter) {
    out <- matrix(NA, nrow=nrow(totals0), ncol=nrow(data0))
    totals0 <- as.data.frame(totals0)
    for (i in 1:nrow(totals0) ) {
      ex <- paste("out[i,] <- as.integer(", sep="")
      for ( k in seq_along(parameter) ) {
        if ( k > 1 ) {
          ex <- paste(ex, "&", sep="")
        }
        ex <- paste(ex, 'data0[,"',parameter[k],'", with=F] =="',totals0[i, parameter[k]],'"', sep="")
      }
      ex <- paste(ex, ")", sep="")
      eval(parse(text=ex))
    }
    out
  }

  inp <- indices_by_constraint(data0, totals0, params$parameter)
  current_totals <- as.numeric(totals0$N)
  weights <- as.integer(data0$weights)

  hh_info <- list()
  hh_info$hh_ids <- as.integer(data0$hid)
  hh_info$hh_head <- as.integer(data0$pid)
  hh_info$hh_head[hh_info$hh_head>1] <- 0L
  hh_info$hh_size <- as.integer(data0$hhsize_calculated)

  w <- .Call("synthPop_calibPop_work", inp=inp, totals=current_totals,
        weights=weights, hh_info=hh_info, params=params, package="synthPop")
  invisible(w)
}


calibPop <- function(data, totals, hid, parameter, split, temp = 30, eps.factor = 0.05, maxiter=200, temp.cooldown = 0.975, factor.cooldown = 0.85, min.temp = 10^-3, verbose=FALSE) {
  x <- NULL
  params <- list()
  params$temp <- as.numeric(temp)[1]
  params$eps_factor = as.numeric(eps.factor)[1]
  params$maxiter = as.integer(maxiter)[1]
  params$temp_cooldown = as.numeric(temp.cooldown)[1]
  params$factor_cooldown = as.numeric(factor.cooldown)[1]
  params$min_temp = as.numeric(min.temp)[1]
  params$verbose <- ifelse(verbose, 1L, 0L)
  params$parameter <- parameter

  parallel <- FALSE
  have_win <- Sys.info()["sysname"] == "Windows"
  nr_cores <- detectCores()
  if ( nr_cores > 2 ) {
    parallel <- TRUE
    nr_cores <- nr_cores-1 # keep one core available
  } else {
    parallel <- FALSE
  }

  # calculate hhsize if not existing
  setkeyv(data, hid)
  res <- data[,.N,by=key(data)]
  setkeyv(res, hid)
  setnames(res, c(hid, "hhsize_calculated"))
  data <- res[data]

  ii <- match(parameter, colnames(data))
  data[,ii] <- data[,lapply(.SD,as.factor),.SDcols=parameter]
  rm(ii)

  ## split problem by "split"-factor
  setkeyv(data,split)
  setkeyv(totals, split)
  split.number <- unique(data[,split,with=F])

  if ( parallel ) {
    # windows
    if ( have_win ) {
      cl <- makePSOCKcluster(nr_cores)
      registerDoParallel(cl)
      final_weights <- foreach(x=1:nrow(split.number), .options.snow=list(preschedule=TRUE)) %dopar% {
        calcFinalWeights(
          data0=data[split.number[x]],
          totals0=totals[which(totals[,split,with=F]==as.character(split.number[x][[split]])),],
          params=params
        )
      }
      stopCluster(cl)
    }
    # linux/mac
    if ( !have_win ) {
      final_weights <- mclapply(1:nrow(split.number), function(x) {
        calcFinalWeights(
          data0=data[split.number[x]],
          totals0=totals[which(totals[,split,with=F]==as.character(split.number[x][[split]])),],
          params=params)
      })
    }
  } else {
    final_weights <- lapply(1:nrow(split.number), function(x) {
      calcFinalWeights(
        data0=data[split.number[x]],
        totals0=totals[which(totals[,split,with=F]==as.character(split.number[x][[split]])),],
        params=params)
    })
  }
  # return dataset with new weights
  data$new.weights <- as.integer(unlist(final_weights))
  data$hhsize_calculated <- NULL
  return(data)
}

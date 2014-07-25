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
    totals0[,parameter] <- apply(totals0[,parameter], 2, myfun)
    for ( x in parameter ) {
      data0[[x]] <- myfun(data0[[x]])
    }
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
  weights <- as.integer(data0[[params$weight]])

  hh_info <- list()
  hh_info$hh_ids <- as.integer(data0[[params$hhid]])
  hh_info$hh_head <- rep(0L, nrow(data0))
  index <- which(sapply(data0[[params$pid]], function(x) { unlist(strsplit(x, "[.]"))[2] } ) == "1")
  hh_info$hh_head[index] <- 1L
  hh_info$hh_size <- as.integer(data0[[params$hhsize]])
  hh_info$median_hhsize <- median(hh_info$hh_size[hh_info$hh_head==1], na.rm=TRUE)

  w <- .Call("synthPop_calibPop_work", inp=inp, totals=current_totals,
    weights=weights, hh_info=hh_info, params=params, package="synthPop")
  invisible(w)
}

calibPop <- function(inp, split, temp = 1, eps.factor = 0.05, maxiter=200,
  temp.cooldown = 0.9, factor.cooldown = 0.85, min.temp = 10^-3, verbose=FALSE) {

  if ( class(inp) != "synthPopObj" ) {
    stop("argument 'inp' must be of class 'synthPopObj'!\n")
  }

  x <- NULL
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

  # generate donors
  data2 <- sampHH(data, sizefactor=2, hid=hid, strata=split, hsize=hhsize)
  data2[[params$weight]] <- 0
  data2 <- data2[, which(!grepl(hid, colnames(data2))), with=FALSE]
  cn <- colnames(data2)
  cn[length(cn)] <- hid
  setnames(data2, cn)
  data2 <- data2[,match(colnames(data), cn), with=FALSE]
  data <- rbind(data, data2)

  parallel <- FALSE
  have_win <- Sys.info()["sysname"] == "Windows"
  nr_cores <- detectCores()
  if ( nr_cores > 2 ) {
    parallel <- TRUE
    nr_cores <- nr_cores-1 # keep one core available
  } else {
    parallel <- FALSE
  }

  ii <- match(parameter, colnames(data))
  data[,ii] <- data[,lapply(.SD,as.factor),.SDcols=parameter]
  rm(ii)

  ## split problem by "split"-factor
  setkeyv(data, split)
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
  data <- data[data$new.weights==1,]
  data[[inp@pop@weight]] <- data$new.weights
  data$new.weights <- NULL
  inp@pop@data <- data
  invisible(inp)
}


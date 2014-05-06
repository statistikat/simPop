specify_pop <- function(data, hhid, hhsize=NULL, pid=NULL, strata=NULL, additional=NULL) {
  if ( !class(hhid)=="character" | length(hhid) != 1 | is.na(match(hhid, colnames(data)))) {
    stop("hhid must be a character defining the variable holding household ids and must be of length 1!\n")
  }
  if ( !is.null(strata) ) {
    if ( !class(strata)=="character" | length(strata) != 1 | is.na(match(strata, colnames(data)))) {
      stop("strata must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
  }

  data <- as.data.table(data)
  setkeyv(data, hhid)

  if ( !is.null(hhsize) ) {
    if ( !class(hhsize)=="character" | length(hhsize) != 1 | is.na(match(hhsize, colnames(data)))) {
      stop("strata must be a character defining the variable holding information on stratas and must be of length 1!\n")
    } 
  } else {
    hhsize <- "hhsize"
    sizes <- data[,.N,by=key(data)]
    setnames(sizes, c(key(data), hhsize))
    data <- data[sizes]
  }
  if ( !is.null(pid) ) {
    if ( !class(pid)=="character" | length(pid) != 1 | is.na(match(pid, colnames(data)))) {
      stop("strata must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
  } else {
    pid <- "pid"
    sizes <- data[,.N,by=key(data)]
    setnames(sizes, c(key(data), hhsize))
    data$pid <- paste(data[,hhid], ".",unlist(sapply(sizes$hhsize, function(x) { seq(1, x) })), sep="")
  }

  vars <- c(hhid, hhsize, pid, strata)
  if ( !is.null(additional) ) {
    if ( any(is.na(match(additional, colnames(data)))) ) {
      stop("at least one additional variable was not found in the input data set!\n")
    } else {
      vars <- c(vars, additional)
    }
  }
  data <- data[,vars,with=F]
  invisible(new("popObj", data=data, hhid=hhid, hhsize=hhsize, pid=pid, strata=strata, additional=additional))
}


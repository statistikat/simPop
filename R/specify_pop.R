specify_sample <- function(data, hhid, strata=NULL) {
  if ( !class(hhid)=="character" | length(hhid) != 1 | is.na(match(hhid, colnames(samp)))) {
    stop("hhid must be a character defining the variable holding household ids and must be of length 1!\n")
  }
  if ( !is.null(strata) ) {
    if ( !class(strata)=="character" | length(strata) != 1 | is.na(match(strata, colnames(samp)))) {
      stop("strata must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
  }

  data <- as.data.table(data)
  setkeyv(data, hhid)
  invisible(new("popObj", data=data, hhid=hhid, strata=strata))
}


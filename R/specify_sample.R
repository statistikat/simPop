specify_sample <- function(data, hhid, weight, strata=NULL) {
  if ( !class(hhid)=="character" | length(hhid) != 1 | is.na(match(hhid, colnames(samp)))) {
    stop("hhid must be a character defining the variable holding household ids and must be of length 1!\n")
  }
  if ( !class(weight)=="character" | length(weight) != 1 | is.na(match(weight, colnames(samp)))) {
    stop("weight must be a character defining the variable holding sampling weights and must be of length 1!\n")
  }
  if ( !is.null(strata) ) {
    if ( !class(strata)=="character" | length(strata) != 1 | is.na(match(strata, colnames(samp)))) {
      stop("strata must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
  }

  data <- as.data.table(data)
  setkeyv(data, hhid)
  invisible(new("sampleObj", data=data, hhid=hhid, weight=weight, strata=strata))
}


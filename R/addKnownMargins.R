addKnownMargins <- function(inp, margins) {
  dataP <- inp@pop@data
  margins <- as.data.frame(margins)
  if ( !class(margins) == "data.frame" ) {
    stop("input argument 'margins' must inherit class 'data.frame'!\n")
  }
  if ( any(duplicated(margins)) ) {
    stop("'margins' must not contain duplicated rows!\n")
  }
  if ( !class(margins[,ncol(margins)]) == "numeric" ) {
    stop("last column of input 'margins' must contain the numbers (must be numeric)!\n")
  } else {
    vals <- margins[,ncol(margins)]
    margins <- margins[,-ncol(margins), drop=FALSE]
  }
  if ( !class(inp) == "synthPopObj" ) {
    stop("input argument 'inp' must be of class 'synthPopObj'!\n")
  }
  if ( !all(colnames(margins) %in% colnames(dataP)) ) {
    stop("all variables existing in input 'margins' must also be existing in the
      synthetic population existing in slot 'pop' of input object 'inp'!\n")
  }

  # order: do all levels exist?
  ind <- match(colnames(margins), colnames(dataP))

  frame <- expand.grid(lapply(ind, function(x) {
    if ( is.factor(dataP[[x]]) ) {
      levels(dataP[[x]])
    } else {
      unique(dataP[[x]])
    }
  }))
  colnames(frame) <- colnames(margins)
  frame <- as.data.table(frame)

  margins$N <- vals
  margins <- as.data.table(margins)
  setkeyv(frame, colnames(frame))
  setkeyv(margins, colnames(margins))
  frame <- merge(frame, margins, all.x=TRUE)
  frame <- frame[,lapply(.SD, as.character)]
  frame$N <- as.numeric(frame$N)
    
  ind <- which(is.na(frame$N))
  if ( length(ind) > 0 ) {
    frame$N[ind] <- 0
  }

  if ( !is.null(inp@table) ) {
    message("Note: currently stored marginals/totals are going to be overwritten!\n")
  }
  inp@table <- frame
  invisible(inp)
}


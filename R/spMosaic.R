spMosaic <- function(x, ...) {
  if ( !class(x) == "spTable" ) { 
    stop("input argument 'x' must be of class 'spTable'!\n")
  } 
  tab <- as.table(x)
  dn <- dimnames(tab)
  dn[[length(dn)]] <- c("Sample", "Population")
  names(dn)[length(dn)] <- "Data"
  dimnames(tab) <- dn
  # define local version of 'cotabplot'
  localCotabplot <- function(x, ..., cond, panel) {
    cotabplot(x, cond="Data", ...)
  }
  # produce plot
  localCotabplot(tab, ...)
}

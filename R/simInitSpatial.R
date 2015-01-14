simInitSpatial <- function(simPopObj, additional, region, tspatial) {
  # simplified version of generateValues_distribution
  generateValues_spatial <- function(dataTable, dataPop, params) {
    if( !nrow(dataTable) ) {
      return(character())
    }
    grid <- expand.grid(dataTable[,params$additional])
    
    # percentages
    perc <- dataTable$freq / sum(dataTable$freq)
    
    # draw and exapand
    sizes <- dataPop[,.N, key(dataPop)]
    sim <- rep(sample(grid[,1], nrow(sizes), prob=perc, replace=TRUE), sizes$N)
    sim
  }
  
  x <- NULL
  dataP <- simPopObj@pop
  dataS <- simPopObj@sample
  data_pop <- dataP@data
  data_sample <- dataS@data
  basic <- simPopObj@basicHHvars

  if ( length(additional) != 1 ) {
    stop("currently exactly one additional spatial variable can be generated!\n")
  }
  if ( length(region) != 1 ) {
    stop("argument 'region' must be a character vector of length 1!\n")
  }  
  if ( additional %in% colnames(data_pop) ) {
    stop("variable specified in argument 'additional' already exists in the population!\n")
  }
  if ( ncol(tspatial) != 3 ) {
    stop("please check input 'tspatial'! It must have exactly 3 columns!\n")
  }
  
  freqs <- tspatial[,ncol(tspatial)]
  if ( !is.numeric(freqs) ) {
    stop("last column of input table must contain numeric values!\n")
  }
  tspatial <- tspatial[,-ncol(tspatial), drop=F]
  
  m <- match(additional, colnames(tspatial))
  if ( is.na(m) ) {
    stop("variable specified in argument 'additional' (",additional,") is not available in input table 'tspatial'!\n")
  }
  add <- tspatial[,m]
  
  m <- match(region, colnames(tspatial))
  if ( is.na(m) ) {
    stop("variable specified in argument 'additional' (",region,") is not available in input table 'tspatial'!\n")
  }
  reg <- tspatial[,m]
  
  # check other variables levels
  m <- match(region, colnames(data_pop))
  if ( is.na(m) ) {
    stop("variable listed in argument 'region' is not available in the synthetic population data of of input 'simPopObj'!\n")
  }
  m <- match(region, colnames(data_sample))
  if ( is.na(m) ) {
    stop("variable listed in argument 'region' is not available in the sample dataset of input 'simPopObj'!\n")
  }  
  
  # generation of our table
  tab <- data.frame(reg, add, freqs)
  colnames(tab) <- c(region, additional, "freq")  
  
  for ( i in 1:2 ) {
    a <- sort(unique(as.character(tab[,i])))
    m <- match(colnames(tab)[i], colnames(data_pop))
    b <- sort(unique(as.character(data_pop[[m]])))    
    if ( any(a!=b) ) {
      stop("We fould a problem in variable ", colnames(tspatial)[i],". Values in input table and synthetic population do not match!\n")
    }
  }

  # list indStrata contains the indices of dataP split by region
  N <- nrow(data_pop)
  indStrata <- split(1:N, data_pop[[region]])
  
  # predictor variables
  predNames <- dataP@hhid  # in spatial case, it can only be the hhid
  
  params <- list()
  params$additional <- additional
  
  values <- lapply(levels(data_sample[[region]]), function(x) {
    generateValues_spatial(
      dataTable=subset(tab, tab[,region]==x),
      dataPop=data_pop[indStrata[[x]], predNames, with=F], params)
  })
  
  ## add new categorical variables to data set and return
  data_pop[[additional]] <- unlist(values)
  
  # check
  simPopObj@pop@data <- data_pop
  invisible(simPopObj)
}
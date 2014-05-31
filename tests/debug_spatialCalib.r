rm(list=ls())
library(devtools)
dev_mode(FALSE)
library(synthPop)

# for testing, source all files in R/*
sapply(list.files("R", full.names=TRUE), source)

# generate district variable
# tmp-testing function to add districts for eusilcS data, 
# each household has the same district
gen_districts <- function(inp) {
  # restrict to hh
  a <- inp[!duplicated(inp$db030),c("db030","db040")]
  spl <- split(a, a$db040)
  regions <- unique(inp$db040)
  
  tmpres <- lapply(1:length(spl), function(x) {
    codes <- paste(x, 1:sample(2:9,1), sep="")
    spl[[x]]$district <- sample(codes, nrow(spl[[x]]), replace=TRUE)
    spl[[x]]
  })
  tmpres <- do.call("rbind", tmpres)
  tmpres <- tmpres[,-c(2)]
  out <- merge(inp, tmpres, by.x=c("db030"), by.y="db030", all.x=TRUE)
  invisible(out)  
}

# simplified version of generateValues_distribution
generateValues_spatial <- function(dataSample, dataPop, params) {
  grid <- params$grid
  additional <- params$additional
  w <- params$w
  
  if( !nrow(dataSample) ) {
    return(character())
  }
  
  # percentages
  perc <- tableWt(dataSample[, c(additional), with=F], dataSample[[w]])
  perc <- perc / sum(perc)
  
  # draw and exapand
  sizes <- dataPop[,.N, key(dataPop)]  
  sim <- rep(sample(grid[,1], nrow(sizes), prob=perc, replace=TRUE), sizes$N)
  sim
}

simInitSpatial <- function(synthPopObj, additional, seed=1) {
  
  x <- NULL
  
  dataP <- synthPopObj@pop
  dataS <- synthPopObj@sample
  data_pop <- dataP@data
  data_sample <- dataS@data
  basic <- synthPopObj@basicHHvars
  
  if ( any(additional %in% colnames(data_pop)) ) {
    stop("variables already exist in the population!\n")
  }  
  if ( length(additional) != 1 ) {
    stop("currently exactly one additional spatial variable can be generated!\n")
  }  
  if ( any(is.na(data_sample[[additional]])) ) {
    stop("there must not be missing values in the sample information about the additional variable!\n")
  }
  
  parallel <- FALSE
  have_win <- Sys.info()["sysname"] == "Windows"
  nr_cores <- detectCores()
  if ( nr_cores > 2 ) {
    parallel <- TRUE
    nr_cores <- nr_cores-1 # keep one core available
  } else {
    parallel <- FALSE
  }
  
  ##### initializations
  if ( !missing(seed) ) {
    set.seed(seed)  # set seed of random number generator
  }
  
  if ( !all(basic %in% colnames(data_pop)) ) {
    stop("undefined variables in the population data -> check your input!\n")
  }
  
  # check sample data against additional parameter
  if ( !all(additional %in% colnames(data_sample)) )  {
    stop("undefined variables in the sample data -> check your input!\n")
  }
  
  # observations with missings are excluded from simulation
  exclude <- getExclude(data_sample)
  if ( length(exclude) > 0 ) {
    data_sample <- data_sample[-exclude,]
  }
  
  # variables are coerced to factors
  data_sample <- checkFactor(data_sample, c(dataS@strata, basic, additional))
  data_pop <- checkFactor(data_pop, c(dataP@strata, basic))
  
  # list indStrata contains the indices of dataP split by strata
  N <- nrow(data_pop)
  indStrata <- split(1:N, data_pop[[dataP@strata]])
  
  ##### simulation
  # predictor variables
  predNames <- dataP@hhid  # in spatial case, it can only be the hhid
  
  params <- list()
  params$grid <- expand.grid(lapply(data_sample[,additional, with=F], levels))
  params$additional <- additional
  params$w <- dataS@weight
  if ( parallel ) {
    # windows
    if ( have_win ) {
      cl <- makePSOCKcluster(nr_cores)
      registerDoParallel(cl)
      values <- foreach(x=levels(data_sample[[dataS@strata]]), .options.snow=list(preschedule=TRUE)) %dopar% {
        generateValues_spatial(
          dataSample=data_sample[data_sample[[dataS@strata]] == x,],
          dataPop=data_pop[indStrata[[x]], predNames, with=F], params
        )
      }
      stopCluster(cl)
    }
    # linux/max
    if ( !have_win ) {
      values <- mclapply(levels(data_sample[[dataS@strata]]), function(x) {
        generateValues_spatial(
          dataSample=data_sample[data_sample[[dataS@strata]] == x,],
          dataPop=data_pop[indStrata[[x]], predNames, with=F], params)
      })
    }
  } else {
    values <- lapply(levels(data_sample[[dataS@strata]]), function(x) {
      generateValues_spatial(
        dataSample=data_sample[data_sample[[dataS@strata]] == x,],
        dataPop=data_pop[indStrata[[x]], predNames, with=F], params)
    })
  }
  
  ## add new categorical variables to data set and return
  data_pop[[additional]] <- unlist(values)

  synthPopObj@pop@data <- data_pop
  invisible(synthPopObj)
}

# todo:
# variable "db040" (region) is available in sample and population
# we have some (distributional, marginal information) about c(vars, district) with district being 
# more detailed spatial information

# add district info
data(eusilcS)
eusilcS <- gen_districts(eusilcS)

inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
synthPopObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))

# Tspatial
additional <- "district"
out <- simInitSpatial(synthPopObj, additional=additional)


par(mfrow=c(1,2))
xx <- table(out@pop@data[[additional]])
yy <- table(out@sample@data[[additional]])
barplot(yy/sum(yy))
barplot(xx/sum(xx))
dev.off()

par(mfrow=c(3,3))
datP <- out@pop@data
datS <- out@sample@data
datP$district <- as.character(datP$district)
splP <- split(datP, datP$db040)
splS <- split(datS, datS$db040)

outP <- lapply(splP, function(x) {
  tab <- table(x[[additional]])
  tab / sum(tab)
})
outS <- lapply(splS, function(x) {
  tab <- table(x[[additional]])
  tab / sum(tab)
})

for ( i in 1:length(outP) ) {
  barplot(cbind(outP[[i]], outS[[i]]), horiz=TRUE)
}
dev.off()



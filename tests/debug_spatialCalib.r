rm(list=ls())
library(devtools)
dev_mode(FALSE)
library(synthPop)
library(data.table)
library(parallel)

# for testing, source all files in R/*
sapply(list.files("R", full.names=TRUE), source)

# generate district variable
# tmp-testing function to add districts for eusilcS data,
# each household has the same district
gen_districts <- function(inp, sample=TRUE) {
  if ( sample ) {
    hhid <- "db030"
    region <- "db040"
  } else {
    hhid <- "hid"
    region <- "region"
  }
  # restrict to hh
  a <- inp[!duplicated(inp[,hhid]),c(hhid, region)]
  spl <- split(a, a[,region])
  regions <- unique(inp[,region])

  tmpres <- lapply(1:length(spl), function(x) {
    codes <- paste(x, 1:sample(2:9,1), sep="")
    spl[[x]]$district <- sample(codes, nrow(spl[[x]]), replace=TRUE)
    spl[[x]]
  })
  tmpres <- do.call("rbind", tmpres)
  tmpres <- tmpres[,-c(2)]
  out <- merge(inp, tmpres, by.x=c(hhid), by.y=hhid, all.x=TRUE)
  invisible(out)
}

simInitSpatial <- function(synthPopObj, additional, region, tspatial, seed=1) {
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
  dataP <- synthPopObj@pop
  dataS <- synthPopObj@sample
  data_pop <- dataP@data
  data_sample <- dataS@data
  basic <- synthPopObj@basicHHvars

  if ( any(additional %in% colnames(data_pop)) ) {
    stop("variables already exists in the population!\n")
  }
  if ( length(additional) != 1 ) {
    stop("currently exactly one additional spatial variable can be generated!\n")
  }

  freqs <- tspatial[,ncol(tspatial)]
  if ( !is.numeric(freqs) ) {
    stop("last column of input table must contain numeric values!\n")
  }
  tspatial <- tspatial[,-ncol(tspatial), drop=F]

  m <- match(additional, colnames(tspatial))
  if ( is.na(m) ) {
    stop("variable specified in argument 'additional' (",additional,") is not available in input table!\n")
  }
  add <- tspatial[,m]
  tspatial <- tspatial[,-m, drop=F]

  # check other variables levels
  m <- match(colnames(tspatial), colnames(data_pop))
  if ( any(is.na(m)) ) {
    stop("some variables listed in the input table are not available in the synthetic population data!\n")
  }

  for ( i in 1:ncol(tspatial) ) {
    a <- sort(unique(as.character(data_pop[[m[i]]])))
    b <- sort(unique(as.character(tspatial[,i])))
    if ( any(a!=b) ) {
      stop("problem in variable ", colnames(tspatial)[i],". Values in input table and synthetic population do not match!\n")
    }
  }

  # generation of our table
  tab <- cbind(tspatial, add, freqs)
  colnames(tab)[(ncol(tab)-1):ncol(tab)] <- c(additional, "freq")

  aggtab <- aggregate(tab$freq, list(tab[,region], tab[,additional]), sum)
  colnames(aggtab) <- c(region, additional, "freq")


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

  # list indStrata contains the indices of dataP split by strata
  N <- nrow(data_pop)
  indStrata <- split(1:N, data_pop[[dataP@strata]])

  # predictor variables
  predNames <- dataP@hhid  # in spatial case, it can only be the hhid

  params <- list()
  params$additional <- additional

  values <- lapply(levels(data_sample[[dataS@strata]]), function(x) {
    generateValues_spatial(
      dataTable=subset(aggtab, aggtab[,region]==x),
      dataPop=data_pop[indStrata[[x]], predNames, with=F], params)
  })

  ## add new categorical variables to data set and return
  data_pop[[additional]] <- unlist(values)

  # check
  synthPopObj@pop@data <- data_pop
  invisible(synthPopObj)
}

quick_sa <- function(synthPopObj, tspatial, additional, region) {
  generate_new_solution <- function(codes, hh_ids, factor) {
    hh_ids_unique <- unique(hh_ids)

    nr_draws <- floor(factor)
    sampled_hh <- sample(hh_ids_unique, nr_draws)
    unique_codes <- unique(codes)

    new_codes <- sample(unique_codes, nr_draws, replace=TRUE)
    xx <- data.frame(hhid=sampled_hh, new_codes)
    yy <- data.frame(hhid=hh_ids, codes=codes)
    zz <- merge(yy, xx, all.x=TRUE)
    ii <- which(!is.na(zz$new_codes))
    zz$codes[ii] <- zz$new_codes[ii]
    return(zz$codes)
  }

  calc_obj <- function(tspatial, pop) {
    cn <- colnames(tspatial)[-ncol(tspatial)]
    setkeyv(pop, cn[1])
    setkeyv(pop, cn)
    res <- pop[,.N, by=cn]
    setkeyv(tspatial, cn[1])
    setkeyv(tspatial, cn)

    m <- merge(tspatial, res)
    return(sum(abs(m$Freq - m$N)))
  }

  calc_factor <- function(obj, hh_head, hh_size) {
    sizes <- hh_size[hh_head]
    return(obj * (median(sizes) / 5))
  }

  temp <- 10
  eps.factor <- 0.05
  maxiter <- 200
  temp_cooldown <- 0.975
  factor_cooldown <- 0.85
  min_temp <- 10^-3
  verb <- TRUE
  cooldown <- 0
  change <- 0
  pop <- synthPopObj@pop@data

  setkeyv(pop, synthPopObj@pop@hhid)
  codes <- codes_orig <- pop[[additional]]
  hh_ids <- pop[[synthPopObj@pop@hhid]]

  cn <- colnames(tspatial)[-ncol(tspatial)]
  setkeyv(pop, cn)
  ts <- as.data.table(tspatial)
  setkeyv(ts, cn)

  # calculate objective value based on initial solution
  obj <- calc_obj(ts, pop)
  best_solution <- list()
  best_solution$obj <- obj
  best_solution$codes <- codes
  eps <- eps.factor*obj

  factor <- calc_factor(obj, hh_head, hh_size)
  if ( obj <= eps ) {
    if ( verb ) {
      cat("We have nothing to do and are already finished!\nValue of objective function: ",obj," | (required precision=",eps,")\n")
    }
    return(codes)
  }

  while ( temp > min_temp ) {
    if ( verb ) {
      cat("current temperature:",temp," (minimal temp=",min_temp,")\n")
    }
    counter <- 1
    while ( counter < maxiter ) {
      setkeyv(pop, synthPopObj@pop@hhid)
      new_codes <- generate_new_solution(codes, hh_ids, factor)
      pop[[additional]] <- new_codes
      spl <- split(pop, pop$db030)

      # calculate new objective value based on new solution
      obj_new <- calc_obj(ts, pop)
      obj_new

      if ( verb ) {
        cat("current value of objective function: ",obj_new," (old value:",obj,")\n")
      }
      if ( obj_new <= eps ) {
        obj <- obj_new;
        codes <- new_codes
        change <- change + 1
        if ( verb ) {
          cat("precision reached!\nValue of objective function: ",obj_new," (required precision=",eps,")\n")
        }
        break;
      }
      # if new solution is better than old one, we accept
      if ( obj_new <= obj ) {
        codes <- new_codes
        obj <- obj_new
        change <- change + 1
        if ( obj_new < best_solution$obj ) {
          cat("storing the best solution found until now!\n")
          best_solution$obj <- obj
          best_solution$codes <- codes
        }
      }
      # if new solution is worse than current, we accept based on temp and value
      if ( obj_new > obj ) {
        if ( runif(1) <= exp(-(obj_new-obj)/temp) ) {
          codes <- new_codes
          obj <- obj_new
          change <- change + 1
        }
      }
      counter <- counter + 1
    }
    # decrease temp and decrease factor accordingly
    # decrease temp by a const fraction (simple method used for testing only)
    temp <- temp_cooldown*temp
    factor <- floor(factor_cooldown*factor)
    if ( factor == 0 ) {
      factor <- 1
    }
    cooldown <- cooldown + 1
    if ( (obj_new <= eps) | (cooldown == 500) ) {
      if ( verb ) {
        cat("Required precision or cooldown state reached!\nValue of objective function: ",obj_new," (required precision=",eps,")\n")
      }
      break;
    }
  }
  invisible(codes)
}


# todo:
# variable "db040" (region) is available in sample and population
# we have some (distributional, marginal information) about c(vars, district) with district being
# more detailed spatial information

# input is a population (census table)
# one variable must be the broader "region", available in pop
# one variable must be the "district",
# additional variables possible but must be available in synthetic population

# add district info
data(eusilcS)
data(eusilcP)
eusilcS <- gen_districts(eusilcS, sample=TRUE)
eusilcP <- gen_districts(eusilcP, sample=FALSE)

inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
synthPopObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))

# age in groups
synthPopObj@pop@data$ageG <- cut(synthPopObj@pop@data$age, c(seq(0, 80, 20),100))
synthPopObj@sample@data$ageG <- cut(synthPopObj@sample@data$age, c(seq(0, 80, 20),100))
eusilcP$ageG <- cut(eusilcP$age, c(seq(0, 80, 20),100))

# table
tspatial <- as.data.frame(xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$ageG + eusilcP$district))
colnames(tspatial) <- c("db040", "rb090", "district", "Freq")
colnames(tspatial) <- c("db040", "ageG", "district", "Freq")
additional <- "district"
region <- "db040"

# first step: quick and fast solution: use distribution of districts within larger region only
synthPopObj <- simInitSpatial(synthPopObj, additional, region, tspatial, seed=1)

# second step: use simulated annealing by region (now just for testing using Burgenland
pop <- synthPopObj@pop@data
pop2 <- subset(pop, pop[[region]]=="Burgenland")
synthPopObj@pop@data <- pop2

tspatial <- subset(tspatial, tspatial[[region]]=="Burgenland")

# applying simulated annlealing to burgenland
tmp <- quick_sa(synthPopObj, tspatial, additional, region)

# improve solution
#synthPopObj <- quick_sa(synthPopObj, tspatial, additional)

# Todo:
# 1) fix sim-annealing and check, if it is really the way it should work
# 2) parallelization of sim-annealing if it works
# 3) if 1) and 2) are fixed, implement slot tspatial into synthPopObj and use it
# 4) write a nice wrapper function around everything




###################
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



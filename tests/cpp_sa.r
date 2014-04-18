rm(list=ls())
library(synthPop)
library(Rcpp)

gen_data <- function() {
  data(eusilcP)
  data(eusilcS)
  
  ## synth pop:
  pop <- eusilcP
  colnames(pop)[3] <- "hhsize"
  ## donor data:
  donors <- sampHH(pop=pop, strata="region", hsize="hhsize")
  
  ## microdata (1) and donor (0)
  pop <- rbind(pop, donors)
  pop$weights <- rep(c(1,0), c(nrow(eusilcP), nrow(donors)))
  pop$new.weights <- pop$weights
  
  index <- pop$weights == 1
  tab <- table(pop$region[index],pop$gender[index], pop$hhsize[index])
  
  ## create target marginals
  totals <- tableWt(eusilcS[,c("db040", "rb090","hsize")], weights=eusilcS$rb050)
  totals <- nrow(pop[index])/sum(totals) * totals
  totals <- ceiling(totals)
  totals <- as.data.table(totals)
  setnames(totals, c("region","gender","hhsize","N"))
  
  data <- data.table(pop, keep.rownames=TRUE)
  totals <- totals
  hid <- "hid"
  parameter <- c("gender","hhsize")
  split <- "region"
  temp <- 10
  eps.factor <-  0.05
  maxiter <- 250
  temp.cooldown <- 0.90
  factor.cooldown <- 0.95
  min.temp <- 10^-2
  sample <- TRUE
  parallel <- FALSE
    
  data0 <- subset(data, data$region=="Vienna")
  totals0 <- subset(totals, totals$region=="Vienna")
  tot_indices <- function(data0, totals0, parameter) {
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
  
  library(Rcpp)
  dat <- list()
  dat$inp <- tot_indices(data0, totals0, parameter)
  dat$totals <- as.numeric(totals0$N)
  dat$weights <- as.integer(data0$weights)
  dat$hh_ids <- as.integer(data0$hid)
  dat$hh_head <- as.integer(data0$pid)
  dat$hh_head[dat$hh_head>1] <- 0L
  dat$hh_size <- as.integer(data0$hhsize)
  dat$verbose <- 1L
  save(dat, file="dat.rdata")
  invisible(dat)
}


if ( !file.exists("dat.rdata") ) {
  gen_data()
}
load("dat.rdata")

sourceCpp("src/calibPop.cpp") 
out <- calib_pop(inp=dat$inp, totals=dat$totals, weights=dat$weights, hh_ids=dat$hh_ids, 
      hh_head=dat$hh_head, hh_size=dat$hh_size, verbose=dat$verbose)
library(data.table)
library(simFrame)
## assumption: this data set is approx 2 x population size
data(eusilcP)
eusilcP <- data.table(eusilcP, keep.rownames=TRUE)
### split into two data sets, one for the synthetic population, one as a donor data set:
#donor <- eusilcP[1:20000,]
### synthetische pop (marginals/constraints not fulfilled):
#eusilcP <- eusilcP[20001:nrow(eusilcP),]
weights <- c(rep(1, 30000),rep(0,nrow(eusilcP)-30000))
eusilcP$weights <- weights
eusilcP$new.weights <- weights

## current marginals:
index <- eusilcP$weights == 1
tab <- table(eusilcP$region[index],eusilcP$gender[index], eusilcP$hsize[index]) # noch eine dritte Var dazugeben
tab
## want to have marginals (artifical constructed, just for testing...)
## usually this is "only" a table of marginals
## we now construct this table just based on eusilcS
## doing like this is the known real marginals/constraints 
library(simPopulation)
data(eusilcS)
totals <- tableWt(eusilcS[,c("db040", "rb090","hsize")], weights=eusilcS$rb050)# noch eine dritte Var dazugeben
totals <- nrow(eusilcP[index])/sum(totals) * totals
## switch regions as in eusilcP
o <- c(1,3,8,2,5,7,4,6,9)
totals <- totals[o,,]

## function for swap
## function for evaluation (already written below)
## function for sa




## we now have/need:
## 1. synthetic population U1 (all with weights==1) that do not respect constraints
## 2. a donor data set, U2 (obs with weights==0)
## 3. constraints/marginals 
## 4. objective function
## aim: synthetic population U1 should fulfil 
## the constraints/marginals --> U*
## swap obs from U2 with U1 by changing weights

### howto:
## simple: start with fraction=1/20 of population size 
## in SA go down with fraction up to 1/Nh, Nh ... number of households


## define objective function - total sums
## use totals from eusilcS

## look at sex - region - hsize

obj <- function(data,totals,weights,parameter){
  
  index <- weights == 1
  synthtotals <- table(data[index,parameter,with=FALSE])
  
  sumtotal <- sum(abs(totals-synthtotals))
  return(sumtotal)
  
}


## function for sample weights
## use difference to marginals
## sample weight euqals positive difference between marginals and synth-data divided by whole deviation
sample.weights <- function(data,totals,weights,parameter, donor=TRUE){
  
index <- weights == 1
synthtotals <- table(data[index,parameter,with=FALSE])
if(donor==FALSE){diff.sum <- synthtotals-totals
                }else{diff.sum <- totals-synthtotals}
            
      ## let categories with too many ind be selected by a vary small chance
      ## use abs-total-diff^-1 as testing option
      pindex <- diff.sum<0
      diff.sum[pindex] <- sum(abs(diff.sum))^-1  
      table.weight <- diff.sum/sum(diff.sum)
      ## get weights for every entry
      
      data[,sweight.new:=table.weight[gender,hsize],by=key(data)]
      
return(data$sweight.new)

}

## resample function
## switches a percentage of weights, given by factor, from 0 to 1
## switches full households always
resample <- function(data, factor){
 
 new.weights <- data$weights
 active <- new.weights == 1  
 inactive <- new.weights == 0 
  
 randomactive <- sample(data$hid[active],factor,prob=data$sweight.synth[active])
 randomactive <- sort(unique(randomactive))
 activeindex <- chmatch(as.character(randomactive),as.character(data$hid))
 activehhsize <- as.numeric(data$hsize[activeindex])
 activeindex <- unlist(sapply(c(1:length(randomactive)), function(x) c(activeindex[x]:c((activeindex+activehhsize-1))[x])))
 
 randominactive <- sample(data$hid[inactive],factor,prob=data$sweight.donor[inactive])  
 randominactive <- sort(unique(randominactive)) 
 inactiveindex <- chmatch(as.character(randominactive),as.character(data$hid)) 
 inactivehhsize <- as.numeric(data$hsize[inactiveindex])
 inactiveindex <- unlist(sapply(c(1:length(randominactive)), function(x) c(inactiveindex[x]:c((inactiveindex+inactivehhsize-1))[x])))
 
 
 new.weights[activeindex] <- 0 
 new.weights[inactiveindex] <- 1
 return(new.weights)
}

  
  
## data...data containing the households ect
## totals...marginals to be reached using simulated annealing
## temp...startin temperatur (use 30 for testing purpose)
## eps...acceptance level for obejctive function (use 500 for testing purpose)
## maxiter...max number of iterations during each temperature setting
## split...problem will be split into sub-problems based of the factorization of split
## parameter...vector which holds the attributs to be considered in the sa-algorithm   

simannealing <- function(data, totals, parameter=c("gender","hsize"), split=c("region"), temp = 30, eps.factor = 0.05, maxiter=200, temp.cooldown = 0.975, factor.cooldown = 0.85, min.temp = 10^-5, sample = FALSE){
  
    
  ## split problem by "split"-factor
  split.number <- unique(data[,which(c(colnames(data) %in% split)==TRUE),with=FALSE])
   
  setkeyv(data,split)
  
  ## initialize for output for every smaller problem    
  solution.list <- c()
  objective.list <- c()
  weights.list <- list()
  s.weights.list <- list()
  factor.list <- c()
  temp.list <- c()
  cooldown.list <- c()
  changes.list <- c()
  
  ## use sa-algorithm for every smaller problem individually
  for(i in 1:nrow(split.number)){
    
    temp0 <- temp
    ##count number of temp coolings - for testing purpose
    cooldown <- 0 
    ## c count if sampleweights are swapped
    change <- 0
    
    ##split data and marginals accordingly
    ##as well as inital subset of weights==1
    data0 <- data[split.number[i]]
      
    ##adjust for different factors
    totals0 <- totals[i,,]
    setkey(data0,"hid")
    
    ## initialize weights-vector in accordance to small problem
    index0 <- data0[round(sum(totals0)),]$hid
    index1 <- max(which(c(data0$hid==index0)==TRUE))
    weights <- c(rep(1, index1),rep(0,nrow(data0)-index1))
    data0$weights <- weights
    data0$new.weights <- weights
    
    ## acceptable error
    eps <- eps.factor*sum(data0$weights)
    ## initial number of swaps per resample
    ## decreases with temperature
    factor <- index1*(1/10)
  

    setkeyv(data0,parameter)
    #parameter <- which(c(colnames(data) %in% split)==TRUE)
    objective <- obj(data = data0 , totals = totals0, weights = data0$weights, parameter = parameter)
    if(sample==TRUE){data0$sweight.synth <- sample.weights(data = data0, totals = totals0 , weights = data0$weights , parameter = parameter, donor=FALSE)
                     data0$sweight.donor <- sample.weights(data = data0, totals = totals0 , weights = data0$weights , parameter = parameter, donor=TRUE)
                      }else{data0$sweight.synth <- c(rep(1,nrow(data0)))
                            data0$sweight.donor <- c(rep(1,nrow(data0)))}
  
    if(objective<=eps){
    
      text <- paste("synthetic-Data fullfils marginals to the Error of", eps)
      solution.list <- c(solution.list,text)
      objective.list <- rbind(objective.list,c(objective,eps))
      weights.list[[i]] <- data0$weights
      s.weights.list[[i]] <- cbind(data0$sweight.synth,data0$sweight.donor)
      factor.list <- c(factor.list,factor)
      temp.list <- c(temp.list,temp0)
      cooldown.list <- c(cooldown.list,cooldown)
      changes.list <- c(changes.list,change)
    }else{
  
   

    ## minimal temperature 
    ## if temperature falls below minimal temp -> terminate
      while(temp0 > min.temp){
      
        n <- 1
            while(n<maxiter){
        
              ## k count if resample does not give a better solution
              #k <- 0

          
              ## swap zeros to one in weights.new
              setkey(data0,"hid")
              data0$new.weights <- resample(data=data0, factor=factor)
              setkeyv(data0,parameter)
              objective.new <- obj(data=data0 , totals=totals0, weights=data0$new.weights,parameter = parameter)
        
              ## if new sample fullfils marginals -> terminate
              if(objective.new <= eps){
                objective <- objective.new
                data0$weights <- data0$new.weights
                change <- change+1
                break
              }
        
              ## choose wether to accepts the resample
              if(objective.new <= objective){ 
                data0$weights <- data0$new.weights 
                objective <- objective.new
                if(sample==TRUE){data0$sweight.synth <- sample.weights(data = data0, totals = totals0 , weights = data0$weights , parameter = parameter, donor=FALSE)
                                 data0$sweight.donor <- sample.weights(data = data0, totals = totals0 , weights = data0$weights , parameter = parameter, donor=TRUE)
                }
                change <- change+1
              }
        
              if(objective.new > objective) {
          
                prob <- exp(-(objective.new-objective)/temp)
                x <- sample(c(0,1), 1,prob=c(1-prob,prob))
          
                if(x == 1){ 
                  data0$weights <- data0$new.weights 
                  objective <- objective.new
                  if(sample==TRUE){data0$sweight.synth <- sample.weights(data = data0, totals = totals0 , weights = data0$weights , parameter = parameter, donor=FALSE)
                                   data0$sweight.donor <- sample.weights(data = data0, totals = totals0 , weights = data0$weights , parameter = parameter, donor=TRUE)
                  }
                  change <- change+1
                }
              
                # count for not changeing weights
                #else {k <- k+1}
          
              }
        
              n <- n+1
        
              ##if temp and factor level do not produce better solutions after 100 iterations cool down temp and reduce factor
              #if(k == 150){break}
            }
      ## decrease temp and decrease factor accordingly
      ## decrease temp by a const fraction (simple method used for testing only)
      temp0 <- temp.cooldown*temp0
      factor <- floor(factor.cooldown*factor)
      if (factor == 0){factor <- 1}
      cooldown <- cooldown + 1
      if(objective.new <= eps | cooldown == 500){break}
       
      }

    


    if(objective <= eps){
      text <- paste("synthetic-Data fullfils marginals to the Error of", eps,"- Synthetic-Data indicated by wheights")
    }else{
      text <- paste("Simulated Annealing failed to optimize objective function to the Error of", eps,"- Synthetic-Data indicated by wheights")  
    }
    
    solution.list <- c(solution.list,text)
    objective.list <- rbind(objective.list,c(objective,eps))
    weights.list[[i]] <- data0$weights
    s.weights.list[[i]] <- cbind(data0$sweight.synth,data0$sweight.donor)
    factor.list <- c(factor.list,factor)
    temp.list <- c(temp.list,temp0)
    cooldown.list <- c(cooldown.list,cooldown)
    changes.list <- c(changes.list,change)
    }
  }
  
  
  solution.list <- cbind(solution.list,split.number)
  objective.list <- cbind(objective.list,split.number)
  factor.list <- cbind(factor.list,split.number)
  temp.list <- cbind(temp.list,split.number)
  cooldown.list <- cbind(cooldown.list,split.number)
  changes.list <- cbind(changes.list,split.number)
  
  
  result <- list(Solution = solution.list , Obj.value = objective.list, sample.weight = s.weights.list, weights = weights.list, parameter = parameter, factor = factor.list, factor.cooldown = factor.cooldown, cooldown = cooldown.list, temp.cooldown = temp.cooldown, Temperature = temp.list , number.changes = changes.list)
  return(result)
  

}


system.time(test0 <- simannealing(data=eusilcP, totals=totals, parameter=c("gender","hsize"), split=c("region"), temp = 10, eps.factor = 0.05, maxiter=250, temp.cooldown = 0.90, factor.cooldown = 0.95, min.temp = 10^-2, sample = TRUE))
test0$Solution
test0$Obj.value
test0$factor
test0$Temperature
test0$cooldown
test0$number.changes
system.time(test1 <- simannealing(data=eusilcP, totals=totals, parameter=c("gender","hsize"), split=c("region"), temp = 10, eps.factor = 0.05, maxiter=250, temp.cooldown = 0.90, factor.cooldown = 0.95, min.temp = 10^-2, sample = FALSE))
test1$Solution
test1$Obj.value
test1$factor
test1$Temperature
test1$cooldown
test1$number.changes
system.time(test2 <- simannealing(data=eusilcP, totals=totals, parameter=c("gender","hsize"), split=c("region"), temp = 10, eps.factor = 0.025, maxiter=250, temp.cooldown = 0.90, factor.cooldown = 0.95, min.temp = 10^-3, sample = TRUE))
test2$Solution
test2$Obj.value
test2$factor
test2$Temperature
test2$cooldown
test2$number.changes
system.time(test3 <- simannealing(data=eusilcP, totals=totals, parameter=c("gender","hsize"), split=c("region"), temp = 10, eps.factor = 0.025, maxiter=250, temp.cooldown = 0.90, factor.cooldown = 0.95, min.temp = 10^-3, sample = FALSE))
test3$Solution
test3$Obj.value
test3$factor
test3$Temperature
test3$cooldown
test3$number.changes
system.time(test2 <- simannealing(data=eusilcP, totals=totals, parameter=c("gender","hsize"), split=c("region"), temp = 10, eps.factor = 0.025, maxiter=250, temp.cooldown = 0.90, factor.cooldown = 0.95, min.temp = 10^-3, sample = TRUE))
test2$Solution
test2$Obj.value
test2$factor
test2$Temperature
test2$cooldown
test2$number.changes
system.time(test3 <- simannealing(data=eusilcP, totals=totals, parameter=c("gender","hsize"), split=c("region"), temp = 10, eps.factor = 0.025, maxiter=250, temp.cooldown = 0.90, factor.cooldown = 0.95, min.temp = 10^-3, sample = FALSE))
test3$Solution
test3$Obj.value
test3$factor
test3$Temperature
test3$cooldown
test3$number.changes















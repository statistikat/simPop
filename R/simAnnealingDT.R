#' @importFrom wrswoR sample_int_crank
#' @importFrom stats sd
#######################################################################################
# Help-functions for simulated annealing without c++
# tries to avoid large memory allocations
dteval <- function(...,envir=parent.frame()){
  eval(parse(text=paste0(...)), envir=envir)
}

# check objective function
checkObjective <- function(totals0,epsH,epsP,epsMinN=0){
  
  namesPos <- grepl("pers",names(totals0))
  hhcond <- perscond <- TRUE
  if(any(!namesPos)){
    hhcond <- sapply(totals0[!namesPos],function(z){
      all(abs(z[["Freq"]]-z[["FreqPop"]])<max(epsMinN,epsH*z[["Freq"]]))
    })
  }
  if(any(namesPos)){
    perscond <- sapply(totals0[!namesPos],function(z){
      all(abs(z[["Freq"]]-z[["FreqPop"]])<max(epsMinN,epsP*z[["Freq"]]))
    })
  }
  
  return(all(hhcond)&all(perscond))
}

# update totals0 with population totals
updateTotals <- function(totals0,data0,hhid,numberPop="weight_choose"){
  
  for(i in 1:length(totals0)){
    
    byVars <- intersect(colnames(totals0[[i]]),colnames(data0))
    if(grepl("pers",names(totals[i]))){
      totals0[[i]] <- merge(totals0[[i]][,mget(c(byVars,"Freq"))],
                 data0[,.(FreqPop=sum(get(numberPop))),by=c(byVars)],
                 by=c(byVars))
    }else{
      totals0[[i]] <- merge(totals0[[i]][,mget(c(byVars,"Freq"))],
                 data0[firstPersonInHousehold==TRUE,.(FreqPop=sum(get(numberPop))),by=c(byVars)],
                 by=c(byVars))
    }
  }
  
  return(totals0)
}

# calc objective function
# absolute sum of number of people/households in population - number of people/household in contingenca table
calcObjective <- function(totals0){
  
  objective <- sapply(totals0,function(z){
    z[,sum(abs(FreqPop-Freq))]
  })
  return(objective)
}

# set number of people to redraw
setRedraw <- function(objective,med_hh=1,fac=2/3){
  
  namesPos <- grepl("pers",names(objective))
  
  if(any(namesPos)){
    val <- max(objective[namesPos])
  }else{
    val <- max(objective[!namesPos])
  }
  
  redraw <- ceiling(val/med_hh*fac)
  return(redraw)
}

# set redraw gap
setRedrawGap <- function(totals0,med_hh=1,scale.redraw=0.5){
  
  med_hh_help <- rep(1,length(totals0))
  med_hh_help[grepl("hh",names(totals0))] <- med_hh
  
  Gap <- rep(0,length(totals0))
  for(i in 1:length(Gap)){
    Gap[i] <- totals0[[i]][,sum(FreqPop-Freq)/med_hh_help[i]]
  }

  redraw_gap <- mean(Gap)*scale.redraw
  return(redraw_gap)  
}


# get probabilites for resampling
getProbabilities <- function(totals0,data0,select_add,select_remove){
  
  # cat("prep totals_diff\n")
  totals_diff <- copy(totals0)
  totals_diff_merged <- NULL
  # cat("start loop\n")
  for(i in seq_along(totals_diff)){
    prob_add <- paste0("prob_add",i)
    prob_remove <- paste0("prob_remove",i)
    # cat(i,"\n")
    totals_diff[[i]][,diff:=as.numeric(Freq-FreqPop)]
    totals_diff[[i]][,c(prob_add):=diff]
    totals_diff[[i]][get(prob_add)<=0,c(prob_add):=exp(sum(get(prob_add)))]
    totals_diff[[i]][,c(prob_remove):=diff*-1]
    totals_diff[[i]][get(prob_remove)<=0,c(prob_remove):=exp(sum(get(prob_remove)))]
    totals_diff[[i]][,helpMergeIndex:=1] # help for merging tables with no common variables
    totals_diff[[i]][,c("Freq","FreqPop","diff"):=NULL]
    if(!is.null( totals_diff_merged)){
      byVar <- intersect(colnames(totals_diff_merged),colnames(totals_diff[[i]]))
      totals_diff_merged <- merge(totals_diff_merged,totals_diff[[i]],by=byVar,allow.cartesian=TRUE,all=TRUE)
    }else{
      totals_diff_merged <- copy(totals_diff[[i]])
    }
  }
  # cat("done\n")
  cnames <- colnames(totals_diff_merged)
  getCols <- cnames[grepl("^prob_add",cnames)]
  totals_diff_merged[,prob_add:=rowMeans(.SD),.SDcols=c(getCols)]
  totals_diff_merged[,c(getCols):=NULL]
  getCols <- cnames[grepl("^prob_remove",cnames)]
  totals_diff_merged[,prob_remove:=rowMeans(.SD),.SDcols=c(getCols)]
  totals_diff_merged[,c(getCols):=NULL]
  
  keyVars <- colnames(totals_diff_merged)[!grepl("prob_remove|prob_add|helpMergeIndex",colnames(totals_diff_merged))]
  
  addIndex <- ((select_add-1)%%nrow(data0)) + 1
  removeIndex <- ((select_remove-1)%%nrow(data0)) + 1

  # cat("merge with data0\n")
  prob_add <- totals_diff_merged[data0[addIndex,mget(keyVars)],prob_add,on=c(keyVars)]
  prob_remove <- totals_diff_merged[data0[removeIndex,mget(keyVars)],prob_remove,on=c(keyVars)]
  probSample <- list(add=prob_add,remove=prob_remove)
  return(probSample)
}


# compare different objectives
compareObjectives <- function(objective,objective_new,med_hh){
  namesPos <- grepl("pers",names(objective))
  
  diffPers <- objective[namesPos]-objective_new[namesPos]
  diffHH <- objective[!namesPos]-objective_new[!namesPos]
  diffAll <- c(diffPers,diffHH)
  w <- c(rep(1,sum(namesPos)),rep(med_hh,sum(!namesPos)))
  
  weighted.mean(diffAll,w) 
}

simAnnealingDT <- function(data0,totals0,params,sizefactor=2,
                           sample.prob=TRUE,choose.temp=FALSE,choose.temp.factor=0.2,
                           scale.redraw=.5,split=NULL,split.level=NULL,observe.times=50,observe.break=0.05){
  N <- V1 <- sim_ID <- weight_choose <- weight_choose_new <- NULL
  ######################################
  ## define variables from param
  
  indTabPers <- which(grepl("pers",names(totals0)))
  indTabHH <- which(grepl("hh",names(totals0)))
  
  epsP <- params[["epsP_factor"]]
  epsH <- params[["epsH_factor"]]
  nd <- nrow(data0)
  hhid <- params[["hhid"]]
  min_temp <- params[["min_temp"]]
  factor_cooldown <- params[["factor_cooldown"]]
  temp_cooldown <- params[["temp_cooldown"]]
  maxiter <- params[["maxiter"]]
  temp <- params[["temp"]]
  hhsize <- params[["hhsize"]]
  parameter <- params[["parameter"]]
  epsMinN <- params[["epsMinN"]]
  npers <- length(indTabPers)
  nhh <- length(indTabHH)
  
  # parameters used for c++ code
  # set index for original order
  dteval("data0[,firstPersonInHousehold:=!duplicated(",hhid,")]")
  data0[,sim_ID:=.I]
  
  ## initialize other parameter
  size_all <- sizefactor+1
  max_n <- size_all* nd
  med_hh <- dteval("data0[!duplicated(",hhid,"),median(as.numeric(as.character(",hhsize,")))]")
  
  if(npers>0){
    init_n <- sapply(totals0[indTabPers],function(z){sum(z[["Freq"]])/med_hh})
  }else{
    init_n <- sapply(totals0[indTabHH],function(z){sum(z[["Freq"]])})
  }
  
  init_n <- mean(init_n)
  
  cooldown <- 0 
  setkeyv(data0,hhid)
  id <- dteval("data0[,",hhid,"]")
  size <- dteval("as.numeric(as.character(data0[,",hhsize,"]))")
  
  
  # choose starting temperatur as percentage of objective function
  if(choose.temp){
    eps <- unlist(c(epsP,epsH))
    eps <- eps[eps>0]
    temp <- max(temp,mean(eps)*choose.temp.factor)
  }
  
  ######################################
  # initialize weights
  init_index <- sample(1:nd,init_n,replace=FALSE)
  
  init_weight <- rep(0L,max_n)
  init_weight[init_index] <- 1
  init_weight <-  updateVecC(init_weight,add_index=init_index, remove_index=c(nd), hhsize=size, hhid=id, sizefactor=size_all)
  data0[ ,weight_choose:=sumVec(init_weight,size_all)]
  
  totals0 <- updateTotals(totals0=totals0,data0=data0,hhid=hhid)
  
  ######################################
  # evaluate objective
  objective <- calcObjective(totals0=totals0)
  
  # define redraw with initial objective value
  redraw <- setRedraw(objective,med_hh=med_hh)
  
  # observe updating of objective function
  # if solution does not improve -> terminate
  # observe only if observe.times>0 and observe.break>0
  observe.count <- 0
  do.observe <- observe.times>0&observe.break>0
  if(do.observe){
    observe.obj <- matrix(0,ncol=npers+nhh,nrow=observe.times)
  }
  
  cat(paste0("Starting simulated Annealing for ",split," ",split.level,"\n"))
  ######################################
  # apply simulated annealing
  if ( checkObjective(totals0,epsH,epsP,epsMinN) ) {
    setkeyv(data0,"sim_ID")
    selectVars <- c(hhid,params[["pid"]],"weight_choose")
    out <- data0[,..selectVars]
    cat(paste0("Convergence successfull for ",split," ",split.level),"\n")
  } else {
    
    ## if objective not fullfilled continue with simannealing
    ## if temperature falls below minimal temp -> terminate
    while( temp > min_temp ) {      
      n <- 1
      while( n<maxiter) {
        # cat("n=",n,"\n")
        # scale redraw for add and remove to keep synthetic totals stable
        
        # cat("set redrawgap")
        redraw_gap <- setRedrawGap(totals0 = totals0,med_hh = med_hh, scale.redraw = scale.redraw)
        # cat("done")
        redraw_add <- max(ceiling(redraw-redraw_gap),1)
        redraw_remove <- max(ceiling(redraw+redraw_gap),1)
        
        # indices to add or remove persons
        select_add <- which(init_weight==0)
        select_remove <- which(init_weight==1)
        
        if(sample.prob){
          #####################################
          # resample
          # get weights for resampling
          # cat("get probabilities\n")
          probs <- getProbabilities(totals0=totals0,data0=data0,select_add=select_add,select_remove=select_remove)
          
          select_add <- select_add[probs[["add"]]>0]
          
          probs[["add"]] <- probs[["add"]][probs[["add"]]>0]
          select_remove <- select_remove[probs[["remove"]]>0]
          probs[["remove"]] <- probs[["remove"]][probs[["remove"]]>0]
          n_add <- length(select_add)
          n_remove <- length(select_remove)
          
          # cat("draw sample\n")
          if(n_add>0){
            add_hh <- select_add[sample_int_crank(n_add,
                                                  min(c(redraw_add,n_add)),
                                                  prob=probs[["add"]])]-1
          }else{
            add_hh <- sample(select_01[[1]],redraw_add)
            
          }
          if(n_remove>0){
            remove_hh <- select_remove[sample_int_crank(n_remove,
                                                        min(c(redraw_remove,n_remove)),
                                                        prob=probs[["remove"]])]-1
          }else{
            remove_hh <- sample(select_01[[2]],redraw_remove)
          }
          
        }else{
          prob_remove <- prob_add <- NULL
          
          ######################################
          # resample
          add_hh <- sample(select_add,n_add)
          remove_hh <- sample(select_remove,n_remove)
        }
        
        ####################################
        ## create new composition
        init_weight_new <- copy(init_weight)

        init_weight_new <-  updateVecC(init_weight_new,add_index=add_hh, remove_index=remove_hh, hhsize=size, hhid=id, sizefactor=size_all)
        data0[ ,weight_choose_new:=sumVec(init_weight_new,size_all)]
        
        ######################################
        # calculate objective
        
        totals0_new <- copy(totals0)
        totals0_new <- updateTotals(totals0=totals0_new,data0=data0,hhid=hhid,numberPop="weight_choose_new")
        objective_new <- calcObjective(totals0_new)
        
        # cat("compare results\n")
        ######################################
        ## if new sample fullfils marginals -> terminate
        if ( checkObjective(totals0,epsH,epsP,epsMinN) ) {
          objective <- objective_new
          data0[,weight_choose:=weight_choose_new]
          break
        }
        
        ######################################
        ## choose wether to accepts the resample
        diffObj <- compareObjectives(objective,objective_new,med_hh=med_hh)
        if ( diffObj>=0 ) { 
          objective <- objective_new
          data0[,weight_choose:=weight_choose_new]
          totals0 <- copy(totals0_new)
          init_weight <- init_weight_new
          
          # update observe variables
          if(do.observe){
            observe.count <- observe.count +1
            if(observe.count>=observe.times){
              
              breakCond <- matrixStats::colSds(observe.obj)/matrixStats::colMeans2(observe.obj)
              breakCond <- max(breakCond)
              if(breakCond< observe.break){
                break # if objective doesnt move anymore break up loop
              }else{
                observe.count <- 1
                observe.obj[observe.count,] <- objective
              }
            }else{
              observe.obj[observe.count,] <- objective
            }
          }
        }
        
        ######################################
        # accept if solution got worse with small probability
        if ( diffObj < 0 ) {      
          prob <- exp(diffObj/temp)
          x <- sample(c(0,1), 1,prob=c(1-prob,prob))
          
          if ( x == 1 ) { 
            objective <- objective_new
            data0[,weight_choose:=weight_choose_new]
            totals0 <- copy(totals0_new)
            init_weight <- init_weight_new
            
            # update observe variables
            if(do.observe){
              observe.count <- observe.count +1
              if(observe.count>=observe.times){
                breakCond <- matrixStats::colSds(observe.obj)/matrixStats::colMeans2(observe.obj)
                breakCond <- max(breakCond)
                if(breakCond< observe.break){
                  break # if objective doesnt move anymore break up loop
                }else{
                  observe.count <- 1
                  observe.obj[observe.count,] <- objective
                }
              }else{
                observe.obj[observe.count,] <- objective
              }
            }
          }
        }    
        n <- n+1
      }
      cat("cooldown\n")
      ## decrease temp and decrease factor accordingly
      ## decrease temp by a const fraction (simple method used for testing only)
      temp <- temp_cooldown*temp
      redraw <- floor(factor_cooldown*redraw)
      if ( redraw == 0 ) {
        redraw <- 1
      }
      cooldown <- cooldown + 1
      if(cooldown%%10==0){
        cat(paste0("Cooldown number ",cooldown,"\n"))
      }
      if ( checkObjective(totals0,epsH,epsP,epsMinN) | cooldown == 500 | redraw<2) {
        break
      }
    }

    # check if convergence was successfull
    if(!checkObjective(totals0,epsH,epsP,epsMinN)){
      cat(paste0("Convergence NOT successfull for ",split," ",split.level),"\n")
    }else{
      cat(paste0("Convergence successfull for ",split," ",split.level),"\n")
    }
    setkeyv(data0,"sim_ID")
    selectVars <- c(hhid,params[["pid"]],"weight_choose")
    out <- data0[,..selectVars] 
  }
  
  out[,c(split):=split.level]
  
  return(out)
}

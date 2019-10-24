#' @importFrom wrswoR sample_int_crank
#' @importFrom stats sd
#######################################################################################
# Help-functions for simulated annealing without c++
# tries to avoid large memory allocations
dteval <- function(...,envir=parent.frame()){
  eval(parse(text=paste0(...)), envir=envir)
}

# check objective function
checkObjective <- function(objective,epsH,epsP){
  
  hhcond <- perscond <- FALSE
  if(!is.null(objective$hh)){
    hhcond <- mapply(`>`,epsH,objective$hh)
  }
  if(!is.null(objective$pers)){
    perscond <- mapply(`>`,epsP,objective$pers)
  }
  return(all(hhcond)&all(perscond))
}

# update totals0 with population totals
updateTotals <- function(totals0,data0,hhid,numberPop="weight_choose"){
  
  if("pers"%in%names(totals0)){
    sapply(totals0$pers,function(z){
      byVars <- intersect(colnames(z),colnames(data0))
      byVars <- byVars[byVars!="unionCode"]
      set(z,j="FreqPop",value=0.0)
      z[data0[,sum(get(numberPop)),by=c(byVars)],FreqPop:=as.numeric(V1),on=c(byVars)]
      return(NULL)
    })
  }
  if("hh"%in%names(totals0)){
    sapply(totals0$hh,function(z){
      byVars <- intersect(colnames(z),colnames(data0))
      byVars <- byVars[byVars!="unionCode"]
      set(z,j="FreqPop",value=0.0)
      z[data0[firstPersonInHousehold==TRUE,sum(get(numberPop)),by=c(byVars)],FreqPop:=as.numeric(V1),on=c(byVars)]
      return(NULL)
    })
  }
  
  return(NULL)
}

# calc objective function
# absolute sum of number of people/households in population - number of people/household in contingenca table
calcObjective <- function(totals0){
  
  hhDiff <- persDiff <- NULL
  if("pers"%in%names(totals0)){
    persDiff <- sapply(totals0$pers,function(z){
      z[,sum(abs(FreqPop-Freq))]
    })
  }
  if("hh"%in%names(totals0)){
    hhDiff <- sapply(totals0$hh,function(z){
      z[,sum(abs(FreqPop-Freq))]
    })
  }
  
  objective <- list(pers=persDiff,hh=hhDiff)
  return(objective)
}

# set number of people to redraw
setRedraw <- function(objective,med_hh=1,fac=2/3){
  if(!is.null(objective$pers)){
    val <- mean(unlist(objective$pers))
  }else{
    val <- mean(unlist(objective$hh))
  }
  
  redraw <- ceiling(val/med_hh*fac)
  return(redraw)
}

# set redraw gap
setRedrawGap <- function(totals0,med_hh=1,scale.redraw=0.5){
  
  hhGap <- persGap <- NULL
  if("pers"%in%names(totals0)){
    persGap <- sapply(totals0$pers,function(z){
      z[,sum(FreqPop-Freq)/med_hh]
    })
  }
  if("hh"%in%names(totals0)){
    hhGap <- sapply(totals0$hh,function(z){
      z[,sum(FreqPop-Freq)/1]
    })
  }

  redraw_gap <- c(hhGap,persGap)
  redraw_gap <- mean(redraw_gap)*scale.redraw
  return(redraw_gap)  
}


# get probabilites for resampling
getProbabilities <- function(totals0,data0,select_add,select_remove){
  totals_diff <- copy(unlist(totals0,recursive = FALSE))
  for(i in seq_along(totals_diff)){
    prob_add <- paste0("prob_add",i)
    prob_remove <- paste0("prob_remove",i)
    
    totals_diff[[i]][,diff:=Freq-FreqPop]
    totals_diff[[i]][,c(prob_add):=diff]
    totals_diff[[i]][get(prob_add)<=0,c(prob_add):=exp(sum(get(prob_add)))]
    totals_diff[[i]][,c(prob_remove):=diff*-1]
    totals_diff[[i]][get(prob_remove)<=0,c(prob_remove):=exp(sum(get(prob_remove)))]
    totals_diff[[i]][,c("Freq","FreqPop","diff"):=NULL]
  }
  totals_diff <- Reduce(function(...) merge(..., all = TRUE,allow.cartesian=TRUE), totals_diff)
  cnames <- colnames(totals_diff)
  getCols <- cnames[grepl("^prob_add",cnames)]
  totals_diff[,prob_add:=matrixStats::rowMaxs(as.matrix(.SD)),.SDcols=c(getCols)]
  totals_diff[,c(getCols):=NULL]
  getCols <- cnames[grepl("^prob_remove",cnames)]
  totals_diff[,prob_remove:=matrixStats::rowMaxs(as.matrix(.SD)),.SDcols=c(getCols)]
  totals_diff[,c(getCols):=NULL]
  
  keyVars <- colnames(totals_diff)[!grepl("prob_remove|prob_add|unionCode",colnames(totals_diff))]
  
  addIndex <- ((select_add-1)%%nrow(data0)) + 1
  prob_add <- totals_diff[data0[addIndex,..keyVars],prob_add,on=c(keyVars)]
  removeIndex <- ((select_remove-1)%%nrow(data0)) + 1
  prob_remove <- totals_diff[data0[removeIndex,..keyVars],prob_add,on=c(keyVars)]
  
  probSample <- list(add=prob_add,remove=prob_remove)
  return(probSample)
}

# compare different objectives
compareObjectives <- function(objective,objective_new,med_hh){
  diffPers <- objective$pers-objective_new$pers
  diffHH <- objective$hh-objective_new$hh
  diffAll <- c(diffPers,diffHH)
  w <- c(rep(1,length(objective$pers)),rep(med_hh,length(objective$hh)))
  
  weighted.mean(diffAll,w) 
}

simAnnealingDT <- function(data0,totals0,params,sizefactor=2,
                           sample.prob=TRUE,choose.temp=FALSE,choose.temp.factor=0.2,
                           scale.redraw=.5,split.level=NULL,observe.times=50,observe.break=0.05){
  N <- V1 <- sim_ID <- weight_choose <- weight_choose_new <- NULL
  ######################################
  ## define variables from param
  
  epsP <- lapply(totals0$pers,function(z){z[,sum(Freq)]*params[["epsP_factor"]]})
  epsH <- lapply(totals0$hh,function(z){z[,sum(Freq)]*params[["epsH_factor"]]})
  nd <- nrow(data0)
  hhid <- params[["hhid"]]
  min_temp <- params[["min_temp"]]
  factor_cooldown <- params[["factor_cooldown"]]
  temp_cooldown <- params[["temp_cooldown"]]
  maxiter <- params[["maxiter"]]
  temp <- params[["temp"]]
  hhsize <- params[["hhsize"]]
  parameter <- params[["parameter"]]
  npers <- length(totals0$pers)
  nhh <- length(totals0$hh)
  
  
  # parameters used for c++ code
  # set index for original order
  dteval("data0[,firstPersonInHousehold:=!duplicated(",hhid,")]")
  data0[,sim_ID:=.I]
  
  ## initialize other parameter
  size_all <- sizefactor+1
  max_n <- size_all* nd
  med_hh <- dteval("data0[!duplicated(",hhid,"),median(as.numeric(as.character(",hhsize,")))]")
  
  if("hh"%in%names(totals0)){
    init_n <- sapply(totals0$hh,function(z){sum(z[["Freq"]])})
  }else{
    init_n <- sapply(totals0$pers,function(z){sum(z[["Freq"]])/med_hh})
  }
  
  init_n <- mean(init_n)

  cooldown <- 0 
  setkeyv(data0,hhid)
  id <- dteval("data0[,",hhid,"]")
  size <- dteval("as.numeric(as.character(data0[,",hhsize,"]))")
  
  # if(sample.prob==TRUE){
  #   init_group <- rep(data0[,ID_GRP],size_all)# used for internal loop
  # }
  
  # choose starting temperatur as percentage of objective function
  if(choose.temp){
    eps <- unlist(c(epsP,epsH))
    eps <- eps[eps>0]
    temp <- max(temp,mean(eps)*choose.temp.factor)
    #min_temp <- temp*temp_cooldown^50
  }
  
  ######################################
  # initialize weights
  init_index <- sample(1:nd,init_n,replace=FALSE)
  init_weight <- rep(0L,max_n)
  init_weight <-  updateVecC(init_weight,add_index=init_index, remove_index=c(max_n), hhsize=size, hhid=id, sizefactor=size_all)
  data0[ ,weight_choose:=sumVec(init_weight,size_all)]
  
  # init_weight <- sample(c(rep(1L,init_n),rep(0L,max_n-init_n)))
  # choose_hh <- matrix(init_weight,nrow=nd,ncol=size_all)
  # adjust choose_hh and select all householdmembers
  # data0[,weight_choose:=matrixStats::rowSums2(choose_hh)]
  # data0[,weight_choose:=max(weight_choose),by=c(hhid)]
  
  updateTotals(totals0=totals0,data0=data0,hhid=hhid)
  
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

  cat(paste0("Starting simulated Annealing for ",split.level,"\n"))
  ######################################
  # apply simulated annealing
  
  if ( checkObjective(objective,epsH,epsP) ) { 
    out <- rowSums(choose_hh)
    cat(paste0("Convergence successfull for ",split.level),"\n")
  } else {
    
    ## if objective not fullfilled continue with simannealing
    ## if temperature falls below minimal temp -> terminate
    while( temp > min_temp ) {      
      n <- 1
      while( n<maxiter ) {
        
        # scale redraw for add and remove to keep synthetic totals stable
        redraw_gap <- setRedrawGap(totals0 = totals0,med_hh = med_hh, scale.redraw = scale.redraw)
        
        redraw_add <- max(ceiling(redraw-redraw_gap),1)
        redraw_remove <- max(ceiling(redraw+redraw_gap),1)
        
        # indices to add or remove persons
        select_add <- which(init_weight==0)
        select_remove <- which(init_weight==1)
        
        if(sample.prob){
          #####################################
          # resample
          # get weights for resampling
          probs <- getProbabilities(totals0=totals0,data0=data0,select_add=select_add,select_remove=select_remove)
         
          select_add <- select_add[probs[["add"]]>0]
          probs[["add"]] <- probs[["add"]][probs[["add"]]>0]
          select_remove <- select_remove[probs[["remove"]]>0]
          probs[["remove"]] <- probs[["remove"]][probs[["remove"]]>0]
          n_add <- length(select_add)
          n_remove <- length(select_remove)
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
        updateTotals(totals0=totals0_new,data0=data0,hhid=hhid,numberPop="weight_choose_new")
        objective_new <- calcObjective(totals0_new)

        ######################################
        ## if new sample fullfils marginals -> terminate
        if ( checkObjective(objective,epsH,epsP) ) {
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
                if(npers>0){
                  observe.obj[observe.count,1:npers] <- objective$pers
                }
                if(nhh>0){
                  observe.obj[observe.count,max(1,npers+1):c(npers+nhh)] <- objective$hh
                }
              }
            }else{
              if(npers>0){
                observe.obj[observe.count,1:npers] <- objective$pers
              }
              if(nhh>0){
                observe.obj[observe.count,max(1,npers+1):c(npers+nhh)] <- objective$hh
              }
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
                  if(npers>0){
                    observe.obj[observe.count,1:npers] <- objective$pers
                  }
                  if(nhh>0){
                    observe.obj[observe.count,max(1,npers+1):c(npers+nhh)] <- objective$hh
                  }
                }
              }else{
                if(npers>0){
                  observe.obj[observe.count,1:npers] <- objective$pers
                }
                if(nhh>0){
                  observe.obj[observe.count,max(1,npers+1):c(npers+nhh)] <- objective$hh
                }
              }
            }
          }
        }    
        n <- n+1
      }
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
      if ( checkObjective(objective,epsH,epsP) | cooldown == 500 | redraw<2) {
        break
      }  
    }
    if(!checkObjective(objective,epsH,epsP)){
      cat(paste0("Convergence NOT successfull for ",split.level),"\n")
    }else{
      cat(paste0("Convergence successfull for ",split.level),"\n")
    }
    setkeyv(data0,"sim_ID")
    out <- data0[,weight_choose]  
  }  
  return(out)
}

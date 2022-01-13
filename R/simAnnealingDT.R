#' @importFrom wrswoR sample_int_crank
#' @importFrom stats sd
#######################################################################################
# Help-functions for simulated annealing without c++
# tries to avoid large memory allocations
dteval <- function(...,envir=parent.frame()){
  eval(parse(text=paste0(...)), envir=envir)
}

# helpfunction to estimate group totals
helpGrouping <- function(x,varSets,varName="weight_choose"){
  firstPersonInHousehold <- NULL
  rbindlist(lapply(varSets,function(z){
    x[,list(FreqPers=sum(get(varName)),FreqHH=sum(get(varName)[firstPersonInHousehold])),by=c(z)]
  }),use.names=TRUE,fill=TRUE,idcol = "GROUP_pop")
  
}

# compare different objectives
compareObjectives <- function(objective,objective_new,med_hh){
  namesPos <- grepl("pers",names(objective))
  
  diffPers <- objective[namesPos]-objective_new[namesPos]
  diffHH <- objective[!namesPos]-objective_new[!namesPos]
  diffHH <- diffHH*med_hh
  diffAll <- c(diffPers,diffHH)
  return(mean(diffAll)) 
}

simAnnealingDT <- function(data0,totals0,params,sizefactor=2,
                           choose.temp=FALSE,choose.temp.factor=0.2,
                           scale.redraw=.5,split=NULL,split.level=NULL,
                           observe.times=50,observe.break=0.05,
                           n.forceCooldown=200){
  N <- V1 <- sim_ID <- weight_choose <- weight_choose_new <- NULL
  FreqType <- factor_med_hh <- Freq <- FreqPopPers <- FreqPopHH <- FreqPers <- 
    Diff <- RowIndex <- MARGININDEX <- GROUP <- NULL
  ######################################
  ## define variables from param
  # totals0 <- totals0[grepl("hh",names(totals0))]

  indTabPers <- grepl("pers",names(totals0))
  indTabHH <- grepl("hh",names(totals0))
  
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
  npers <- sum(indTabPers)
  nhh <- sum(indTabHH)
  verbose <- params[["verbose"]]
  
  # parameters used for c++ code
  # set index for original order
  dteval("data0[,firstPersonInHousehold:=!duplicated(",hhid,")]")
  # data0[,firstPersonInHousehold:=TRUE]
  data0[,sim_ID:=.I]
  
  # initialize other parameter
  size_all <- sizefactor+1
  max_n <- size_all* nd
  med_hh <- dteval("data0[!duplicated(",hhid,"),median(as.numeric(as.character(",hhsize,")))]")
  
  cooldown <- 0 
  setkeyv(data0,hhid)
  id <- dteval("data0[,",hhid,"]")
  size <- dteval("as.numeric(as.character(data0[,",hhsize,"]))")
  # size <- rep(1,length(size))
  
  # choose starting temperatur as percentage of objective function
  if(choose.temp){
    eps <- sapply(totals0,`[[`,"Freq")
    eps <- c(unlist(eps[indTabPers])*epsP,
             unlist(eps[indTabHH])*epsH)
    eps <- eps[eps>0]
    temp <- max(temp,mean(eps)*choose.temp.factor)
  }
  
  ######################################
  # initialize weights
  data0[,weight_choose:=0]
  data0[list(split.level),weight_choose:=1,on=c(split)]
  init_index <- data0[weight_choose==1,which=TRUE]
  init_weight <- rep(0L,max_n)
  init_weight[init_index] <- 1L
  indexAddRemove <- splitVector(init_weight)

  marginTable <- copy(totals0)
  marginSets <- list()
  for(i in seq_along(marginTable)){
    
    marginSets[[i]] <- colnames(marginTable[[i]])[colnames(marginTable[[i]])!="Freq"]
    if(indTabHH[i]){
      newName <- "FreqHH"
      factor_i <- med_hh
      eps_i <- epsH
    }else{
      newName <- "FreqPers"
      factor_i <- 1
      eps_i <- epsP
    }
    marginTable[[i]][,FreqType:=newName]
    marginTable[[i]][,factor_med_hh:=factor_i]
    marginTable[[i]][,eps:=pmax(eps_i*Freq,epsMinN)]
    
  }
  names(marginSets) <- names(totals0)
  marginNames <- unique(unlist(marginSets))
  marginTable <- rbindlist(marginTable,use.names=TRUE,fill=TRUE,idcol = "GROUP")
  setcolorder(marginTable,marginNames)
  setkeyv(marginTable,marginNames)
  

  popGrouped <- helpGrouping(x=data0[weight_choose>0],varSets=marginSets)
  # marginTable[popGrouped,c("FreqPers","FreqHH"):=list(FreqPopPers,FreqPopHH),on=c(marginNames)]
  marginTable <- merge(marginTable,popGrouped,by=c(marginNames),all=TRUE)
  if(any(is.na(marginTable[["GROUP"]]))){
    # fill up missings if combination exists in population but not in margins
    fill_vars <- c("GROUP", "Freq", "FreqType", "factor_med_hh", "eps")
    marginTable[is.na(GROUP),c("GROUP","Freq","FreqType","eps"):=.(GROUP_pop,
                                                            0,
                                                            fifelse(grepl("pers",GROUP_pop),"FreqPers","FreqHH"),
                                                            epsMinN)]
    marginTable[,factor_med_hh:=factor_med_hh[!is.na(factor_med_hh)][1],by=.(GROUP)]
  }
  marginTable[,GROUP_pop:=NULL]
  marginTable[is.na(FreqPers),c("FreqPers","FreqHH"):=0]
  marginTable[,Diff:=(Freq-get(unlist(.BY))),by=list(FreqType)]
  # marginTable[abs(Diff)<eps,Diff:=sign(Diff)*sqrt(abs(Diff))]
  # marginTable[,Diffprob:=Diff/factor_med_hh]
  marginTable[,RowIndex:=.I-1]
  h_margins <- as.integer(marginTable$FreqType=="FreqHH") # for c++ function needed
  
  ######################################
  # initialize unique groups in data0
  # for efficient calculation of sampling probabilities
  data0Unique <- data0[,list(Npop=.N),by=c(marginNames)]
  data0Unique[,MARGININDEX:=.I-1]
  
  for(i in 1:length(marginSets)){
    newName <- paste0(names(marginSets)[i],"RowIndex")
    margin.i <- marginTable[GROUP==names(marginSets)[i]]
    data0Unique[margin.i,c(newName):=RowIndex,on=c(marginSets[[i]])]
  }
  rowIndices <- colnames(data0Unique)[grepl("RowIndex",colnames(data0Unique))]
  
  indexMatrix <- as.matrix(data0Unique[,rowIndices, with = FALSE])
  indexdata0 <- data0Unique[data0,MARGININDEX,on=c(marginNames)]
  
  ######################################
  # evaluate objective
  objective <- marginTable[,sum(abs(Diff)),by=list(GROUP)][,sqrt(mean(V1^2))]    # *factor_med_hh
  
  # define redraw with initial objective value
  redraw <- marginTable[,max(abs(Diff)*2/3)]
  redrawMax <- sapply(totals0,function(z){
    sum(z[["Freq"]])
  })
  redrawMax[indTabPers] <- redrawMax[indTabPers]/med_hh
  redrawMax <- round(max(redrawMax*0.05))
  redraw <- min(redraw,redrawMax)
  # cat("redraw=",redraw,"\n")
  
  # observe updating of objective function
  # if solution does not improve -> terminate
  # observe only if observe.times>0 and observe.break>0
  observe.count <- 0
  do.observe <- observe.times>0&observe.break>0
  if(do.observe){
    observe.obj <- matrix(0,ncol=npers+nhh,nrow=observe.times)
  }
  
  noChange <- 0
  updatepSet <- TRUE
  number_iterations <- 0
  ######################################
  # apply simulated annealing
  if ( marginTable[,sum(eps)>=sum(abs(Diff)),by=.(GROUP)][,all(V1==TRUE)] ) { # marginTable[,sum(eps)>=sum(abs(Diff))]
    
    if(verbose){
      message(paste0(split," ",split.level," already fulfills error margins\n"))
    }
    
    setkeyv(data0,"sim_ID")
    selectVars <- c(hhid,params[["pid"]],"weight_choose")
    out <- data0[, selectVars, with = FALSE]
  } else {
    
    if(verbose){
      message(paste0("Starting simulated Annealing for ",split," ",split.level,"\n"))
    }
    ## if objective not fullfilled continue with simannealing
    ## if temperature falls below minimal temp -> terminate
    while( temp > min_temp ) {      
      n <- 1
      # seedX <- sample(1:10000,maxiter)
      
      while( n<maxiter) {
        number_iterations <- number_iterations +1
        # scale redraw for add and remove to keep synthetic totals stable
        
        # message("set redrawgap")
        # redraw_gap <- setRedrawGap(totals0 = totals0,med_hh = med_hh, scale.redraw = scale.redraw)
        redraw_gap <- marginTable[,mean.default(Diff)]*scale.redraw
        # message("done")
        redraw_add <- max(ceiling(redraw+redraw_gap),1)
        redraw_remove <- max(ceiling(redraw-redraw_gap),1)
        
        #####################################
        # resample
        # get weights for resampling
        if(updatepSet){
          pSet <- calcProbabilities(indexMat=indexMatrix,x=marginTable[["Diff"]],Npop=data0Unique[["Npop"]],
                                    indexData = indexdata0,initWeight = init_weight,
                                    indexAdd = indexAddRemove[["indexAdd"]],
                                    indexRemove = indexAddRemove[["indexRemove"]],
                                    n_add = redraw_add, n_remove= redraw_remove)
          updatepSet <- FALSE
        }
       
        
        # set.seed(seedX[n])
        if(pSet[["nAdd"]]>0){
          add_hh <- pSet[["indexAdd"]][sample_int_crank(pSet[["nAdd"]],
                                                        min(c(redraw_add,pSet[["nAdd"]])),
                                                        prob=pSet[["probAdd"]])]
        }else{
          add_hh <- sample(indexAddRemove[["indexAdd"]],min(c(redraw_add,length(indexAddRemove[["indexAdd"]]))))
        }

        if(pSet[["nRemove"]]>0){
          remove_hh <- pSet[["indexRemove"]][sample_int_crank(pSet[["nRemove"]],
                                                              min(c(redraw_remove,pSet[["nRemove"]])),
                                                              prob=pSet[["probRemove"]])]
        }else{
          remove_hh <- sample(indexAddRemove[["indexRemove"]],min(c(redraw_add,length(indexAddRemove[["indexRemove"]]))))
        }

        ####################################
        ## create new composition
        # init_weight_new2 <- copy(init_weight)
        # # update init_weight_new for add_hh and remove_hh
        # # and update diff vector from marginTable to account for changes due
        # # to add_hh and remove_hh
        # init_weight_new2 <-  updateVecC(init_weight_new2,add_index=add_hh, remove_index=remove_hh, hhsize=size, hhid=id, sizefactor=size_all)
        # set(data0,i=NULL,j="weight_choose_new",value=sumVec(init_weight_new2,size_all))
        # # ######################################
        # # # calculate new margins 
        # marginTable_new <- copy(marginTable)
        # marginTable_new[,c("FreqPers","FreqHH"):=NULL]
        # popGrouped <- helpGrouping(x=data0[weight_choose_new>0],varSets=marginSets,varName = "weight_choose_new")
        # marginTable_new[popGrouped,c("FreqPers","FreqHH"):=list(FreqPopPers,FreqPopHH),on=c(marginNames)]
        # marginTable_new[is.na(FreqPers),c("FreqPers","FreqHH"):=0]
        # marginTable_new[,Diff:=(Freq-get(unlist(.BY))),by=list(FreqType)]
        
        init_weight_new <- copy(init_weight)
        new_solution <-  updateObjectiveC(init_weight_new,add_index=add_hh, remove_index=remove_hh,
                                          hhsize=size, hhid=id, sizefactor=size_all,
                                          indexMat = indexMatrix, indexData = indexdata0,
                                          diff = copy(marginTable[["Diff"]]),
                                          householdMargin = h_margins)
        # if(any(abs(new_solution[["diff_new"]]-marginTable_new$Diff)>1e-8)){
        #   stop()
        # }
        set(marginTable,j="Diff_new",value=new_solution[["diff_new"]])
        objective_new <- marginTable[,sum(abs(Diff_new)),by=list(GROUP)][,sqrt(mean(V1^2))]  #  *factor_med_hh
        
        if(verbose){
          cat("objective old = ",objective,"  -  objective new = ",objective_new,"\n")
        }
        ######################################
        ## if new sample fullfils marginals -> terminate
        if ( marginTable[,sum(eps)>=sum(abs(Diff_new)),by=.(GROUP)][,all(V1==TRUE)] ) { # marginTable[,all(eps>=abs(Diff))]
          objective <- objective_new
          marginTable[,c("Diff"):=NULL]
          setnames(marginTable,"Diff_new","Diff")
          set(data0,i=NULL,j="weight_choose",value=sumVec(new_solution[["init_weight"]],size_all))
          break
        }
        # data0[weight_choose_new>weight_choose,.N,by=c(marginNames)][order(N)]
        # message("number households to have:",mean(sapply(totals0,function(z){sum(z$Freq)})),"\n")
        # message("number households sampled:",data0[!duplicated(hid),sum(weight_choose_new)],"\n")
        ######################################
        ## choose wether to accepts the resample
        diffObj <- objective - objective_new
        # diffObj
        if ( diffObj>=0 ) {
          if(verbose){
            message("solution improved!\n")
            cat("objective new:",objective_new,"\n\n")
          }
          marginTable[,c("Diff"):=NULL]
          setnames(marginTable,"Diff_new","Diff")
          set(data0,i=NULL,j="weight_choose",value=sumVec(new_solution[["init_weight"]],size_all))
          init_weight <- copy(new_solution[["init_weight"]])
          indexAddRemove <- splitVector(init_weight)
          objective <- objective_new
          updatepSet <- TRUE
          # update observe variables
          if(do.observe){
            observe.count <- observe.count +1
            if(observe.count>=observe.times){
              breakCond <- colSds(observe.obj)/colMeans2(observe.obj)
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
            if(verbose){
              message("accept worse!\n")
              cat("objective new:",objective_new,"\n\n")
            }
            # message("new solution accepted!\n")
            marginTable[,c("Diff"):=NULL]
            setnames(marginTable,"Diff_new","Diff")
            set(data0,i=NULL,j="weight_choose",value=sumVec(new_solution[["init_weight"]],size_all))
            init_weight <- copy(new_solution[["init_weight"]])
            indexAddRemove <- splitVector(init_weight)
            objective <- objective_new
            updatepSet <- TRUE
            # update observe variables
            if(do.observe){
              observe.count <- observe.count +1
              if(observe.count>=observe.times){
                breakCond <- colSds(observe.obj)/colMeans2(observe.obj)
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
        
        # break inner loop if 
        if(updatepSet){
          noChange <- 0
        }else{
          noChange <- noChange + 1
        }
        if(noChange>n.forceCooldown){
          break # break if solution does not move
        }
      }
      ## decrease temp and decrease factor accordingly
      ## decrease temp by a const fraction (simple method used for testing only)
      noChange <- 0
      updatepSet <- TRUE
      temp <- temp_cooldown*temp
      redraw <- floor(factor_cooldown*redraw)
      if ( redraw == 0 ) {
        redraw <- 1
      }
      cooldown <- cooldown + 1
      if(cooldown%%10==0 & verbose){
        message(paste0("Cooldown number ",cooldown,"\n"))
      }
      if ( marginTable[,sum(eps)>=sum(abs(Diff)),by=.(GROUP)][,all(V1==TRUE)] | cooldown == 500 |(cooldown > 50 & redraw<2)) { # marginTable[,all(eps>=abs(Diff))]
        break
      }
    }

    # check if convergence was successfull
    if(verbose){
      if(marginTable[,sum(eps)>=sum(abs(Diff)),by=.(GROUP)][,all(V1==TRUE)]){ # marginTable[,all(eps>=abs(Diff))]
        message(paste0("Convergence successfull for ",split," ",split.level)," in ",number_iterations,"iterations!\n")
      }else{
        message(paste0("Convergence NOT successfull for ",split," ",split.level)," in ",number_iterations," iterations!\n")
      }
    }

    setkeyv(data0,"sim_ID")
    selectVars <- c(hhid,params[["pid"]],
                    params[["redist.var"]],params[["hhid_orig"]],params[["pid_orig"]],  # <- only not NULL if variable was "redistributed"
                    "weight_choose")
    out <- data0[, selectVars, with = FALSE]
  }
  
  out[,c(split):=split.level]

  return(out)
}

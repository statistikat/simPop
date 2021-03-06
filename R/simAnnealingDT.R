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
    x[,list(FreqPopPers=sum(get(varName)),FreqPopHH=sum(get(varName)[firstPersonInHousehold])),by=c(z)]
  }),use.names=TRUE,fill=TRUE)
  
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
                           n.forceCooldown=10){
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
  init_weight[init_index] <- 1
  
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
  marginTable[popGrouped,c("FreqPers","FreqHH"):=list(FreqPopPers,FreqPopHH),on=c(marginNames)]
  marginTable[is.na(FreqPers),c("FreqPers","FreqHH"):=0]
  marginTable[,Diff:=(Freq-get(unlist(.BY))),by=list(FreqType)]
  # marginTable[abs(Diff)<eps,Diff:=sign(Diff)*sqrt(abs(Diff))]
  # marginTable[,Diffprob:=Diff/factor_med_hh]
  marginTable[,RowIndex:=.I-1]
  
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
  objective <- marginTable[eps<abs(Diff),sum(abs(Diff)^2),by=list(GROUP)][,sqrt(mean(V1))]
  
  # define redraw with initial objective value
  redraw <- marginTable[,max(abs(Diff)*2/3)]
  redrawMax <- sapply(totals0,function(z){
    sum(z[["Freq"]])
  })
  redrawMax[indTabPers] <- redrawMax[indTabPers]/med_hh
  redrawMax <- round(max(redrawMax*0.025))
  redraw <- min(redraw,redrawMax)
  
  
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
  message(paste0("Starting simulated Annealing for ",split," ",split.level,"\n"))
  ######################################
  # apply simulated annealing
  if ( marginTable[,all(eps>=abs(Diff))] ) {
    setkeyv(data0,"sim_ID")
    selectVars <- c(hhid,params[["pid"]],"weight_choose")
    out <- data0[, selectVars, with = FALSE]
    message(paste0("Convergence successfull for ",split," ",split.level),"\n")
  } else {
    
    ## if objective not fullfilled continue with simannealing
    ## if temperature falls below minimal temp -> terminate
    while( temp > min_temp ) {      
      n <- 1
      # seedX <- sample(1:10000,maxiter)
      
      while( n<maxiter) {
        
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
        # probAddC <- apply(indexMatrix,1,function(z){
        #   calcCase(x[z+1])
        # })
        # if(any(abs(probAddC-probAdd)>1e-10)){
        #   stop()
        # }
        # 
        # probAdd <- probAdd/data0Unique[["Npop"]]
        # probRemove <- probAdd*-1
        # probAdd[probAdd<=0] <- exp(sum(probAdd[probAdd<=0]))
        # probRemove[probRemove<=0] <- exp(sum(probRemove[probRemove<=0]))
        
        if(updatepSet){
          pSet <- calcProbabilities(indexMat=indexMatrix,x=marginTable[["Diff"]],Npop=data0Unique[["Npop"]],
                                    indexData = indexdata0,initWeight = init_weight,
                                    indexAdd = indexAddRemove[["indexAdd"]],
                                    indexRemove = indexAddRemove[["indexRemove"]])
          updatepSet <- FALSE
        }
 
          
        # data0Unique[,probAdd0:=probAdd]
        # data0Unique[,probRemove0:=probRemove]
        # probAdd <- probAdd[indexdata0+1]
        # probAdd <- probAdd[(indexAddRemove[["indexAdd"]]%%nrow(data0))+1]
        # indexAdd <- indexAddRemove[["indexAdd"]][probAdd>0]
        # probAdd <- probAdd[probAdd>0]
        # 
        # probAdd2 <- pSet$probAdd
        # indexAdd2 <- pSet$indexAdd
        # 
        # probRemove <- probRemove[indexdata0+1]
        # probRemove <- probRemove[(indexAddRemove[["indexRemove"]]%%nrow(data0))+1]
        # indexRemove <- indexAddRemove[["indexRemove"]][probRemove>0]
        # probRemove <- probRemove[probRemove>0]
        # 
        # probRemove2 <- pSet$probRemove
        # indexRemove2 <- pSet$indexRemove
        # 
        # if(any(abs(probAdd2-probAdd)>1e-10)){
        #   stop()
        # }
        # if(any(abs(probRemove2-probRemove)>1e-10)){
        #   stop()
        # }
        # if(!all.equal(indexAdd2,indexAdd)){
        #   stop()
        # }
        # if(!all.equal(indexRemove2,indexRemove)){
        #   stop()
        # }
        # all(abs(probRemove-pSet$probRemove)<1e-10)
        # pSet$nAdd <- 10
        # pSet$nRemove <- 10
        # indexdata0[head(pSet$indexAdd[order(pSet$probAdd,decreasing = TRUE)])]
        # message("draw sample\n")
        
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

        # set.seed(seedX[n])
        # add_hh2 <- indexAdd[sample_int_crank(length(indexAdd),
        #                                     min(c(redraw_add,length(indexAdd))),
        #                                     prob=probAdd)]
        # # probAdd[which(indexAdd%in%add_hh)]
        # remove_hh2 <- indexRemove[sample_int_crank(length(indexRemove),
        #                                     min(c(redraw_remove,length(indexRemove))),
        #                                     prob=probRemove)]
        # if(!all.equal( add_hh2, add_hh)){
        #   stop()
        # } 
        # if(!all.equal(remove_hh2,remove_hh)){
        #   stop()
        # }
        ####################################
        ## create new composition
        init_weight_new <- copy(init_weight)
        init_weight_new <-  updateVecC(init_weight_new,add_index=add_hh, remove_index=remove_hh, hhsize=size, hhid=id, sizefactor=size_all)
        set(data0,i=NULL,j="weight_choose_new",value=sumVec(init_weight_new,size_all))
        
        
        ######################################
        # calculate new margins 
        marginTable_new <- copy(marginTable)
        marginTable_new[,c("FreqPers","FreqHH"):=NULL]
        popGrouped <- helpGrouping(x=data0[weight_choose_new>0],varSets=marginSets,varName = "weight_choose_new")
        marginTable_new[popGrouped,c("FreqPers","FreqHH"):=list(FreqPopPers,FreqPopHH),on=c(marginNames)]
        marginTable_new[is.na(FreqPers),c("FreqPers","FreqHH"):=0]
        marginTable_new[,Diff:=(Freq-get(unlist(.BY))),by=list(FreqType)]
        # marginTable_new[,Diffprob:=Diff/factor_med_hh]
        # marginTable_new[abs(Diff)<eps,Diff:=sign(Diff)*sqrt(abs(Diff))]
        
        objective_new <- marginTable_new[eps<abs(Diff),sum(abs(Diff)^2),by=list(GROUP)][,sqrt(mean(V1))]
        # objective_new
        # objective
        # marginTable_new[,sum(abs(Diff)^2),by=list(GROUP)][,sqrt(mean(V1))]
        # marginTable[,sum(abs(Diff)^2),by=list(GROUP)]
        # message("compare results\n")
        ######################################
        ## if new sample fullfils marginals -> terminate
        if ( marginTable_new[,all(eps>=abs(Diff))] ) {
          objective <- objective_new
          marginTable <-  copy(marginTable_new)
          data0[,weight_choose:=weight_choose_new]
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
          # message("solution improved!\n")
          marginTable <- copy(marginTable_new)
          data0[,weight_choose:=weight_choose_new]
          init_weight <- copy(init_weight_new)
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
            
            # message("new solution accepted!\n")
            marginTable <- copy(marginTable_new)
            data0[,weight_choose:=weight_choose_new]
            init_weight <- copy(init_weight_new)
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
      if(cooldown%%10==0){
        message(paste0("Cooldown number ",cooldown,"\n"))
      }
      if ( marginTable[,all(eps>=abs(Diff))] | cooldown == 500 | redraw<2) {
        break
      }
    }

    # check if convergence was successfull
    if(!marginTable[,all(eps>=abs(Diff))]){
      message(paste0("Convergence NOT successfull for ",split," ",split.level),"\n")
    }else{
      message(paste0("Convergence successfull for ",split," ",split.level),"\n")
    }
    setkeyv(data0,"sim_ID")
    selectVars <- c(hhid,params[["pid"]],"weight_choose")
    out <- data0[, selectVars, with = FALSE]
  }
  
  out[,c(split):=split.level]
  
  return(out)
}

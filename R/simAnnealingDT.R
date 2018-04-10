#' @importFrom wrswoR sample_int_crank
#' @importFrom stats sd
#######################################################################################
# Help-functions for simulated annealing without c++
# tries to avoid large memory allocations
dteval <- function(...,envir=parent.frame()){
  eval(parse(text=paste0(...)), envir=envir)
}

simAnnealingDT <- function(data0,totals0,params,sizefactor=2,
                           sample.prob=TRUE,choose.temp=FALSE,choose.temp.factor=0.2,
                           scale.redraw=.5,split.level=NULL,observe.times=50,observe.break=0.05){
  ID_GRP <- N <- V1 <- sim_ID <- weight_choose <- weight_choose_new <- NULL
  ######################################
  ## define variables from param
  eps <- params[["eps_factor"]]*totals0[,sum(N)]
  nd <- nrow(data0)
  hhid <- params[["hhid"]]
  min_temp <- params[["min_temp"]]
  factor_cooldown <- params[["factor_cooldown"]]
  temp_cooldown <- params[["temp_cooldown"]]
  maxiter <- params[["maxiter"]]
  eps_factor <- params[["eps_factor"]]
  temp <- params[["temp"]]
  hhsize <- params[["hhsize"]]
  parameter <- params[["parameter"]]
  
  # parameters used for c++ code
  # set index for original order
  data0[,sim_ID:=.I]
  
  ## initialize other parameter
  size_all <- sizefactor+1
  max_n <- size_all* nd
  med_hh <- dteval("data0[!duplicated(",hhid,"),median(as.numeric(as.character(",hhsize,")))]")
  init_n <- round(totals0[,sum(N)/med_hh])
  redraw <- ceiling(med_hh/5 *init_n)
  cooldown <- 0 
  totals0[,ID_GRP:=.GRP,by=c(parameter)]
  data0 <- merge(data0,unique(subset(totals0,select=c("ID_GRP",parameter))),by=c(parameter),all.x=TRUE)
  setkeyv(data0,hhid)
  id <- dteval("data0[,",hhid,"]")
  size <- dteval("as.numeric(as.character(data0[,",hhsize,"]))")
  
  if(sample.prob==TRUE){
    init_group <- rep(data0[,ID_GRP],size_all)# used for internal loop
  }
  
  # choose starting temperatur as percentage of objective function
  if(choose.temp){
    temp <- max(temp,eps*choose.temp.factor)
    #min_temp <- temp*temp_cooldown^50
  }
  
  ######################################
  # initialize weights
  init_weight <- c(rep(1L,nd),rep(0L,max_n-nd))
  # init_weight <- sample(c(rep(1L,init_n),rep(0L,max_n-init_n)))
  choose_hh <- matrix(init_weight,nrow=nd,ncol=size_all)
  
  # adjust choose_hh and select all householdmembers
  dteval("choose_hh[data0[,",hhid,"%in%data0[choose_hh[,",1:size_all,"]>0,",hhid,"]],",1:size_all,"] <- 1L")
  data0[,weight_choose:=rowSums(choose_hh)]
  
  ######################################
  # evaluate objective
  totals_diff <- merge(totals0,data0[,sum(weight_choose),by=ID_GRP],by="ID_GRP")
  totals_diff[is.na(V1),V1:=NA]
  totals_diff[,diff:=V1-N]
  
  objective <- totals_diff[,sum(abs(diff))]
  
  # define redraw with initial objective value
  redraw <- ceiling(objective/med_hh*2/3)
  
  # observe updating of objective function
  # if solution does not improve -> terminate
  # observe only if observe.times>0 and observe.break>0
  observe.count <- 1
  do.observe <- observe.times>0&observe.break>0
  if(do.observe){
    observe.obj <- rep(objective,observe.times)
  }

  
  cat(paste0("Starting simulated Annealing for ",split.level,"\n"))
  ######################################
  # apply simulated annealing
  if ( objective <= eps ) { 
    out <- rowSums(choose_hh)
    cat(paste0("Convergence successfull for ",split.level),"\n")
  } else {
    
    ## if objective not fullfilled continue with simannealing
    ## if temperature falls below minimal temp -> terminate
    while( temp > min_temp ) {      
      n <- 1
      while( n<maxiter ) {
        
        # scale redraw for add and remove to keep synthetic totals stable
        redraw_gap <- totals_diff[,c(sum(V1)-sum(N))/med_hh]*scale.redraw
        
        #if(abs(redraw_gap)<redraw){
        redraw_add <- max(ceiling(redraw-redraw_gap),1)
        redraw_remove <- max(ceiling(redraw+redraw_gap),1)
        #}else{
        #  redraw_add <- redraw_remove <- redraw
        #}
        
        
        if(sample.prob){
          # get weights for resampling
          totals_diff[,prob_add:=diff*-1]
          totals_diff[prob_add<=0,prob_add:=exp(sum(prob_add))]
          totals_diff[,prob_remove:=diff]
          totals_diff[prob_remove<=0,prob_remove:=exp(sum(prob_remove))]
          
          #sample_add <- sample(totals_diff[,ID_GRP],redraw,prob=totals_diff[,prob_add],replace=TRUE)
          #sample_remove <- sample(totals_diff[,ID_GRP],redraw,prob=totals_diff[,prob_remove],replace=TRUE)
          
          #####################################
          # resample
          select_01 <- select_equal(x=init_weight,val1=0,val2=1)
          select_add <- select_01[[1]]+1
          select_remove <- select_01[[2]]+1
          prob_add <- totals_diff[list(init_group[select_add]),prob_add]
          prob_remove <- totals_diff[list(init_group[select_remove]),prob_remove]
          
          select_add <- select_add[prob_add>0]
          select_remove <- select_remove[prob_remove>0]
          n_add <- length(select_add)
          n_remove <- length(select_remove)
          if(n_add>0){
            add_hh <- select_add[sample_int_crank(n_add,
                                                  min(c(redraw_add,n_add)),
                                                  prob=prob_add[prob_add>0])]-1
          }else{
            add_hh <- sample(select_01[[1]],redraw_add)
            
          }
          if(n_remove>0){
            remove_hh <- select_remove[sample_int_crank(n_remove,
                                                        min(c(redraw_remove,n_remove)),
                                                        prob=prob_remove[prob_remove>0])]-1
          }else{
            remove_hh <- sample(select_01[[2]],redraw_remove)
          }
          #add_hh <- sample(select_add,redraw_add)-1
          #remove_hh <- sample(select_remove,redraw_remove)-1
        }else{
          prob_remove <- prob_add <- NULL
          
          ######################################
          # resample
          select_01 <- select_equal(x=init_weight,val1=0,val2=1)
          
          add_hh <- sample(select_01[[1]],redraw_add,prob=prob_add)
          remove_hh <- sample(select_01[[2]],redraw_remove,prob=prob_remove)
        }
        
        ####################################
        ## create new composition
        init_weight_new <- copy(init_weight)
        init_weight_new <-  updateVecC(init_weight_new,add_index=add_hh, remove_index=remove_hh, hhsize=size, hhid=id, sizefactor=size_all)
        
        data0[ ,weight_choose_new:=sumVec(init_weight_new,size_all)]
        ######################################
        # calculate objective
        totals_diff_new <- merge(totals0,data0[,sum(weight_choose_new),by="ID_GRP"],by="ID_GRP",all.x=TRUE)
        totals_diff_new[is.na(V1),V1:=0]
        totals_diff_new[,diff:=V1-N]
        objective.new <- totals_diff_new[,sum(abs(diff))]
        
        ######################################
        ## if new sample fullfils marginals -> terminate
        if ( objective.new <= eps ) {
          objective <- objective.new
          data0[,weight_choose:=weight_choose_new]
          totals_diff <- copy(totals_diff_new)
          init_weight <- init_weight_new
          break
        }
        
        ######################################
        ## choose wether to accepts the resample
        if ( objective.new <= objective ) { 
          objective <- objective.new
          data0[,weight_choose:=weight_choose_new]
          totals_diff <- copy(totals_diff_new)
          init_weight <- init_weight_new
          
          # update observe variables
          if(do.observe){
            observe.count <- observe.count +1
            if(observe.count>observe.times){
              if(sd(observe.obj)/mean(observe.obj)< observe.break){
                break # if objective doesnt move anymore break up loop
              }else{
                observe.count <- 1
                observe.obj[observe.count] <- objective
              }
            }else{
              observe.obj[observe.count] <- objective
            }
          }
        }
        
        ######################################
        # accept if solution got worse with small probability
        if ( objective.new > objective ) {      
          prob <- exp(-(objective.new-objective)/temp)
          x <- sample(c(0,1), 1,prob=c(1-prob,prob))
          
          if ( x == 1 ) { 
            objective <- objective.new
            data0[,weight_choose:=weight_choose_new]
            totals_diff <- copy(totals_diff_new)
            init_weight <- init_weight_new
            
            # update observe variables
            if(do.observe){
              observe.count <- observe.count +1
              if(observe.count>observe.times){
                if(sd(observe.obj)/mean(observe.obj)< observe.break){
                  break # if objective doesnt move anymore break up loop
                }else{
                  observe.count <- 1
                  observe.obj[observe.count] <- objective
                }
              }else{
                observe.obj[observe.count] <- objective
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
      if ( objective.new <= eps | cooldown == 500 | redraw<2) {
        break
      }  
    }
    if(objective>eps){
      cat(paste0("Convergence NOT successfull for ",split.level),"\n")
    }else{
      cat(paste0("Convergence successfull for ",split.level),"\n")
    }
    setkeyv(data0,"sim_ID")
    out <- data0[,weight_choose]  
  }  
  return(out)
}

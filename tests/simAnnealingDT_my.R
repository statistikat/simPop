#######################################################################################
# Help-functions for simulated annealing without c++
# tries to avoid large memory allocations
dteval <- function(...,envir=parent.frame()){
  eval(parse(text=paste0(...)), envir=envir)
}

simAnnealingDT <- function(data0,totals0,params,sizefactor=2,sample.prob=TRUE,choose.temp=FALSE){
  
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
  
  # parameters used für c++ code
  # set index for original order
  data0[,sim_ID:=.I]
  setkeyv(data0,hhid)
  size <- dteval("as.numeric(as.character(data0[,",hhsize,"]))")
  id <- dteval("data0[,",hhid,"]")
  
  ## initialize other parameter
  size_all <- sizefactor+1
  max_n <- size_all* nd
  med_hh <- dteval("data0[!duplicated(",hhid,"),median(as.numeric(as.character(",hhsize,")))]")
  init_n <- round(totals0[,sum(N)/med_hh])
  cooldown <- 0 
  totals0[,ID_GRP:=.GRP,by=c(parameter)]
  data0[,ID_GRP:=.GRP,by=c(parameter)]
  if(sample.prob==TRUE){
    init_group <- rep(data0[,ID_GRP],size_all)# used for internal loop
  }
  
  # choose starting temperatur as percentage of objective function
  if(choose.temp){
    temp <- max(temp,eps)
    #min_temp <- temp*temp_cooldown^50
  }
  
  ######################################
  # initialize weights
  # init_weight <- sample(c(rep(1L,init_n),rep(0L,max_n-init_n)))
  init_weight <- c(c(rep(1L,nd),rep(0L,max_n-nd)))
  choose_hh <- matrix(init_weight,nrow=nd,ncol=size_all)
  
  # adjust choose_hh and select all householdmembers
  dteval("choose_hh[data0[,",hhid,"%in%data0[choose_hh[,",1:size_all,"]>0,",hhid,"]],",1:size_all,"] <- 1L")
  data0[,weight_choose:=rowSums(choose_hh)]
  
  ######################################
  # evaluate objective
  totals_diff <- merge(data0[,sum(weight_choose),by=ID_GRP],totals0,by="ID_GRP")
  totals_diff[,diff:=V1-N]
  
  objective <- totals_diff[,sum(abs(diff))]
  
  # define redraw with initial objective value
  redraw <- ceiling(objective*.1)
  draw_bound <- eps/nrow(totals0)
  draw_bound <- 20
  
  ######################################
  # apply simulated annealing
  if ( objective <= eps ) { 
    out <- rowSums(choose_hh) 
  } else {
    
    ## if objective not fullfilled continue with simannealing
    ## if temperature falls below minimal temp -> terminate
    while( temp > min_temp ) {      
      n <- 1
      while( n<maxiter ) {
        if(sample.prob){
          # get weights for resampling
          totals_diff[,prob_add:=diff*-1]
          totals_diff[prob_add<=0,prob_add:=exp(sum(prob_add))]
          totals_diff[,prob_remove:=diff]
          totals_diff[prob_remove<=0,prob_remove:=exp(sum(prob_remove))]
          
          # sample_add <- sample(totals_diff[,ID_GRP],redraw,prob=totals_diff[,prob_add],replace=TRUE)
          # sample_remove <- sample(totals_diff[,ID_GRP],redraw,prob=totals_diff[,prob_remove],replace=TRUE)
          
          # num_add <- tableC(sample_add)
          # group_add <- as.numeric(names(sample_add))
          # num_remove <- tableC(sample_remove)
          # group_remove <- as.numeric(names(sample_remove))
          
          # select_01 <- select_equal(x=init_weight,val1=0,val2=1)
          # select_add <- select_01[[1]]
          # select_remove <- select_01[[2]]
          
          # add_hh <- sample_group(select_add,group_x=init_group[select_add+1],group=group_add,group_num=num_add,replace=FALSE)
          
          #####################################
          # resample
          select_01 <- select_equal(x=init_weight,val1=0,val2=1)
          select_add <- select_01[[1]]+1
          select_remove <- select_01[[2]]+1
          prob_add <- totals_diff[init_group[select_add],prob_add]
          prob_remove <- totals_diff[init_group[select_remove],prob_remove]
          
          select_add <- select_add[prob_add>draw_bound]
          select_remove <- select_remove[prob_remove>draw_bound]
          # select_add <- head(select_add,20)
          # prob_add <- round(head(prob_add,20))
          # prob_add[sample_int_expj(length(select_add),5,prob=prob_add/sum(prob_add))]
          # 
          # select_add <- sample(51700:51799,20)
          # prob_add <- c(sample(80:150,10),rep(0,10))
          # prob_add[sample_int_expj(length(select_add),size=5,prob=prob_add)]
          
          #add_hh <- select_add[sample_int_expj(length(select_add),redraw,prob=prob_add[prob_add>0])]-1
          #remove_hh <- select_remove[sample_int_expj(length(select_remove),redraw,prob=prob_remove[prob_remove>0])]-1
          add_hh <- sample(select_add,redraw)
          remove_hh <- sample(select_remove,redraw)
          
        }else{
          prob_remove <- prob_add <- NULL
          
          ######################################
          # resample
          select_01 <- select_equal(x=init_weight,val1=0,val2=1)
          
          add_hh <- sample(select_01[[1]],redraw,prob=prob_add)
          remove_hh <- sample(select_01[[2]],redraw,prob=prob_remove)
        }
        
        ####################################
        ## create new composition
        init_weight_new <- copy(init_weight)
        init_weight_new <-  updateVecC(init_weight_new,add_index=add_hh, remove_index=remove_hh, hhsize=size, hhid=id, sizefactor=size_all)
        
        data0[ ,weight_choose_new:=sumVec(init_weight_new,size_all)]
        ######################################
        # calculate objective
        totals_diff_new <- merge(data0[,sum(weight_choose_new),by=c(parameter)],totals0,by=parameter)
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
      if ( objective.new <= eps | cooldown == 500 ) {
        break
      }  
    }
    setkeyv(data0,"sim_ID")
    out <- data0[,weight_choose]  
  }  
  return(out)
}

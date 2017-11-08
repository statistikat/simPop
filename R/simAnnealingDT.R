#######################################################################################
# Help-functions for simulated annealing without c++
# tries to avoid large memory allocations
dteval <- function(...,envir=parent.frame()){
  eval(parse(text=paste0(...)), envir=envir)
}

simAnnealingDT <- function(data0,totals0,params,sizefactor=2,sample.prob=FALSE,choose.temp=TRUE){
  
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
  
  ## initialize other parameter
  size_all <- sizefactor+1
  max_n <- size_all* nd
  med_hh <- dteval("data0[!duplicated(",hhid,"),median(as.numeric(as.character(",hhsize,")))]")
  init_n <- round(totals0[,sum(N)/med_hh])
  choose_hh <- matrix(0L,nrow= nd,ncol=size_all)
  redraw <- ceiling(med_hh/5 *init_n)
  cooldown <- 0 
  # choose starting temperatur as percentage of objective function
  if(choose.temp){
    temp <- max(temp,eps)
    min_temp <- temp*temp_cooldown^50
  }
  
  ######################################
  # initialize weights
  init_weight <- sample(c(rep(1L,init_n),rep(0L,max_n-init_n)))
  init_weight <- matrix(init_weight,nrow=nd,ncol=size_all)
  
  # adjust choose_hh and select all householdmembers
  dteval("choose_hh[data0[,",hhid,"%in%data0[init_weight[,",1:size_all,"]>0,",hhid,"]],",1:size_all,"] <- 1L")
  data0[,weight_choose:=rowSums(choose_hh)]
  
  ######################################
  # evaluate objective
  totals_diff <- merge(data0[,sum(weight_choose),by=c(parameter)],totals0,by=parameter)
  totals_diff[,diff:=V1-N]
  
  objective <- totals_diff[,sum(abs(diff))]
  
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
        
        init_weight <- as.vector(choose_hh)
        if(sample.prob){
          # get weights for resampling
          prob_sample <- rep(totals_diff[data0[,mget(parameter)],diff,on=c(parameter)],size_all)
          total_remove <- totals_diff[,diff]
          total_add <- total_remove*-1
          
          # sample households with probabilites derived by differences to marginals
          # should speed up convergence
          neg_add <- sum(abs(total_add[total_add<0]))^-1
          neg_remove <- sum(abs(total_remove[total_remove<0]))^-1
          total_add[total_add<0] <- neg_add
          total_remove[total_remove<0] <- neg_remove
          
          prob_remove <- prob_sample[init_weight==1]
          prob_remove[prob_remove<0] <- neg_remove
          prob_add <- prob_sample[init_weight==0]*-1
          prob_add[prob_add<0] <- neg_add
          
          prob_remove <- prob_remove/sum(total_remove)
          prob_add <- prob_add/sum(total_add)
          
          if(all(prob_add==0)){
            prob_add <- NULL
          }
          if(all(prob_remove==0)){
            prob_remove <- NULL
          }
        }else{
          prob_remove <- prob_add <- NULL
        }

        ######################################
        # resample

        add_hh <- sample(which(init_weight==0),redraw,prob=prob_add)

        remove_hh <- sample(which(init_weight==1),redraw,prob=prob_remove)
        
        init_weight <- matrix(init_weight,nrow=nd,ncol=size_all)
        ######################################
        # remove households
        remove_col <- floor((remove_hh-1)/nd)+1
        remove_row <- remove_hh%%nd
        remove_row[remove_row==0] <- nd
        remove_col_e <- unique(remove_col)
        remove_row_e <- unlist(lapply(remove_col_e,function(z){paste(remove_row[remove_col==z],collapse=",")}))
        remove_row <- unique(remove_row)
        dteval("init_weight[data0[,",hhid,"%in%data0[c(",remove_row_e,"),",hhid,"]],",remove_col_e,"] <- 0")
        # add households
        add_col <- floor((add_hh-1)/nd)+1
        add_row <- add_hh%%nd
        add_row[add_row==0] <- nd
        add_col_e <- unique(add_col)
        add_row_e <- unlist(lapply(add_col_e,function(z){paste(add_row[add_col==z],collapse=",")}))
        add_row <- unique(add_row)
        dteval("init_weight[data0[,",hhid,"%in%data0[c(",add_row_e,"),",hhid,"]],",add_col_e,"] <- 1")
        
        data0[ ,weight_choose_new:=rowSums(init_weight)]
        ######################################
        # calculate objective
        totals_diff <- merge(data0[,sum(weight_choose_new),by=c(parameter)],totals0,by=parameter)
        totals_diff[,diff:=V1-N]
        totals_diff[,sum(abs(diff))]
        objective.new <- totals_diff[,sum(abs(diff))]
        
        ######################################
        ## if new sample fullfils marginals -> terminate
        if ( objective.new <= eps ) {
          objective <- objective.new
          data0[ ,weight_choose:=weight_choose_new]
          break
        }
        
        ######################################
        ## choose wether to accepts the resample
        if ( objective.new <= objective ) { 
          objective <- objective.new
          data0[ ,weight_choose:=weight_choose_new]
        }
        
        ######################################
        # accept if solution got worse with small probability
        if ( objective.new > objective ) {      
          prob <- exp(-(objective.new-objective)/temp)
          x <- sample(c(0,1), 1,prob=c(1-prob,prob))
          
          if ( x == 1 ) { 
            objective <- objective.new
            data0[ ,weight_choose:=weight_choose_new]
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
    
    out <- data0[,weight_choose]  
  }  
  return(out)
}
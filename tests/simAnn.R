# set.seed(1234)
# test <- rep(FALSE,500)
# 
# for(t in 1:500){
# 
#   init_weight <- as.vector(choose_hh)
#   
#   add_hh <- sample(which(init_weight==0),redraw,prob=prob_add)
#   
#   remove_hh <- sample(which(init_weight==1),redraw,prob=prob_remove)
#   
#   init_weight <- matrix(init_weight,nrow=nd,ncol=size_all)
#   init_weight_old <- init_weight_C <- init_weight
#   ######################################
#   # remove households
#   remove_col <- floor(remove_hh/(nd+1))
#   remove_row <- remove_hh%%nd-1
#   remove_row[remove_row==-1] <- nd-1
#   
#   add_col <- floor(add_hh/(nd+1))
#   add_row <- add_hh%%nd-1
#   add_row[add_row==-1] <- nd-1
#   
#   init_weight_new <- updateMatC(M=init_weight_C,add_row=add_row,add_col=add_col,
#                                 remove_row=remove_row,remove_col=remove_col,hhsize=size,hhid=id)
#   
#   ######################################
#   # remove households
#   remove_col <- floor(remove_hh/(nd+1))+1
#   remove_row <- remove_hh%%nd
#   remove_row[remove_row==0] <- nd
#   remove_col_e <- unique(remove_col)
#   remove_row_e <- unlist(lapply(remove_col_e,function(z){paste(remove_row[remove_col==z],collapse=",")}))
#   remove_row <- unique(remove_row)
#   dteval("init_weight_old[data0[",hhid,"%in%data0[c(",remove_row_e,"),",hhid,"],which=TRUE],",remove_col_e,"] <- 0")
#   
#   # add households
#   add_col <- floor(add_hh/(nd+1))+1
#   add_row <- add_hh%%nd
#   add_row[add_row==0] <- nd
#   add_col_e <- unique(add_col)
#   add_row_e <- unlist(lapply(add_col_e,function(z){paste(add_row[add_col==z],collapse=",")}))
#   add_row <- unique(add_row)
#   dteval("init_weight_old[data0[",hhid,"%in%data0[c(",add_row_e,"),",hhid,"],which=TRUE],",add_col_e,"] <- 1")
#   
#   if(!all(init_weight_old==init_weight_new)){
#     break
#   }else{
#     test[t] <- all(init_weight_old==init_weight_new)
#   }
# }
# all(test==TRUE)
# 
# funC <- function(remove_hh,add_hh,init_weight,nd,size,id){
#   ######################################
#   # remove households
#   remove_col <- floor(remove_hh/(nd+1))
#   remove_row <- remove_hh%%nd-1
#   remove_row[remove_row==-1] <- nd-1
#   
#   add_col <- floor(add_hh/(nd+1))
#   add_row <- add_hh%%nd-1
#   add_row[add_row==-1] <- nd-1
#   
#   init_weight_new <- updateMatC(M=init_weight,add_row=add_row,add_col=add_col,
#                                 remove_row=remove_row,remove_col=remove_col,hhsize=size,hhid=id)
# }
# 
# 
# funDT <- function(remove_hh,add_hh,init_weight_old,nd,hhid){
#   ######################################
#   # remove households
#   remove_col <- floor(remove_hh/(nd+1))+1
#   remove_row <- remove_hh%%nd
#   remove_row[remove_row==0] <- nd
#   remove_col_e <- unique(remove_col)
#   remove_row_e <- unlist(lapply(remove_col_e,function(z){paste(remove_row[remove_col==z],collapse=",")}))
#   remove_row <- unique(remove_row)
#   dteval("init_weight_old[data0[",hhid,"%in%data0[c(",remove_row_e,"),",hhid,"],which=TRUE],",remove_col_e,"] <- 0")
#   
#   # add households
#   add_col <- floor(add_hh/(nd+1))+1
#   add_row <- add_hh%%nd
#   add_row[add_row==0] <- nd
#   add_col_e <- unique(add_col)
#   add_row_e <- unlist(lapply(add_col_e,function(z){paste(add_row[add_col==z],collapse=",")}))
#   add_row <- unique(add_row)
#   dteval("init_weight_old[data0[",hhid,"%in%data0[c(",add_row_e,"),",hhid,"],which=TRUE],",add_col_e,"] <- 1")
#   
# }

# microbenchmark(funDT(remove_hh,add_hh,init_weight_old,nd,hhid),
#                funC(remove_hh,add_hh,init_weight,nd,size,id))
# 
# 
# # create dummy example
# nd <- 50
# id <- sort(sample(1:15,nd,replace=TRUE))
# hsize <- rep(0,length(id))
# for(i in 1:15){
#   hsize[i==id] <- sum(i==id)
# }
# 
# data_new <- data.table(hhid=id,hsize)
# 
# init_weight_t <- matrix(sample(c(0,1),length(id)*4,replace = TRUE,prob=c(.9,.1)),nrow=length(id),ncol=4)
# 
# for(i in 1:4){
#   setid <- unique(id[init_weight_t[,i]==1])
#   init_weight_t[id%in%setid,i] <- 1
# }
# 
# set.seed(1234)
# test <- rep(FALSE,10000)
# for(t in 1:10000){
#   add_hh <- sample(which(init_weight_t==0),5,prob=prob_add)
#   
#   remove_hh <- sample(which(init_weight_t==1),5,prob=prob_remove)
#   
#   init_weight_C <- matrix(init_weight_t,nrow=nd,ncol=4)
#   init_weight_DT <- init_weight_C
#   ######################################
#   # remove households
#   remove_col <- floor(remove_hh/(nd+1))
#   remove_row <- remove_hh%%nd-1
#   remove_row[remove_row==-1] <- nd-1
#   
#   add_col <- floor(add_hh/(nd+1))
#   add_row <- add_hh%%nd-1
#   add_row[add_row==-1] <- nd-1
#   
#   init_weight_C<- updateMatC(M=init_weight_C,add_row=add_row,add_col=add_col,
#                              remove_row=remove_row,remove_col=remove_col,hhsize=hsize,hhid=id)
#   
#   
#   ######################################
#   # remove households
#   remove_col <- floor(remove_hh/(nd+1))+1
#   remove_row <- remove_hh%%nd
#   remove_row[remove_row==0] <- nd
#   remove_col_e <- unique(remove_col)
#   remove_row_e <- unlist(lapply(remove_col_e,function(z){paste(remove_row[remove_col==z],collapse=",")}))
#   remove_row <- unique(remove_row)
#   dteval("init_weight_DT[data_new[",hhid,"%in%data_new[c(",remove_row_e,"),",hhid,"],which=TRUE],",remove_col_e,"] <- 0")
#   
#   # add households
#   add_col <- floor(add_hh/(nd+1))+1
#   add_row <- add_hh%%nd
#   add_row[add_row==0] <- nd
#   add_col_e <- unique(add_col)
#   add_row_e <- unlist(lapply(add_col_e,function(z){paste(add_row[add_col==z],collapse=",")}))
#   add_row <- unique(add_row)
#   dteval("init_weight_DT[data_new[",hhid,"%in%data_new[c(",add_row_e,"),",hhid,"],which=TRUE],",add_col_e,"] <- 1")
#   
#   test[t] <- all(init_weight_C==init_weight_DT)
# }
# 
# all(test==TRUE)
# get probabilites for resampling
getProbabilities <- function(totals0,data0,select_add,select_remove,med_hh){
  helpMergeIndex <- Freq <- FreqPop <- NULL
  # cat("prep totals_diff\n")
  totals_diff <- copy(totals0)
  totals_diff_merged <- NULL
  # cat("start loop\n")
  for(i in seq_along(totals_diff)){
    prob_add <- paste0("prob_add",i)
    prob_remove <- paste0("prob_remove",i)
    # cat(i,"\n")
    set(totals_diff[[i]], j = "diff", value =
          as.numeric(totals_diff[[i]][["Freq"]] - totals_diff[[i]][["FreqPop"]]
          ))
    if(grepl("hh",names(totals_diff)[i])){
      set(totals_diff[[i]], j = "diff", value =
            totals_diff[[i]][["diff"]]*med_hh)
    }
    # help for merging tables with no common variables
    set(totals_diff[[i]], j = "helpMergeIndex", value = 1)
    set(totals_diff[[i]], j = prob_add, value = totals_diff[[i]][["diff"]])
    set(totals_diff[[i]], j = prob_remove, value = -totals_diff[[i]][["diff"]])
    set(totals_diff[[i]], j = c("Freq","FreqPop","diff"), value = NULL)
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
  totals_diff_merged[,prob_add:=rowSums(.SD),.SDcols=c(getCols)]
  totals_diff_merged[prob_add<=0,prob_add:=exp(sum(prob_add))]
  set(totals_diff_merged, j = getCols, value = NULL)
  getCols <- cnames[grepl("^prob_remove",cnames)]
  totals_diff_merged[,prob_remove:=rowSums(.SD),.SDcols=c(getCols)]
  totals_diff_merged[prob_remove<=0,prob_remove:=exp(sum(prob_remove))]
  set(totals_diff_merged, j = getCols, value = NULL)
  
  keyVars <- colnames(totals_diff_merged)[!grepl("prob_remove|prob_add|helpMergeIndex",colnames(totals_diff_merged))]
  
  addIndex <- ((select_add-1)%%nrow(data0)) + 1
  removeIndex <- ((select_remove-1)%%nrow(data0)) + 1
  
  # cat("merge with data0\n")
  prob_add <- totals_diff_merged[data0[addIndex,mget(keyVars)],prob_add,on=c(keyVars)]
  prob_remove <- totals_diff_merged[data0[removeIndex,mget(keyVars)],prob_remove,on=c(keyVars)]
  probSample <- list(add=prob_add,remove=prob_remove)
  return(probSample)
}

# update totals0 with population totals
updateTotals <- function(totals0,data0,hhid,numberPop="weight_choose"){
  firstPersonInHousehold <- NULL
  for(i in 1:length(totals0)){
    
    byVars <- intersect(colnames(totals0[[i]]),colnames(data0))
    if(grepl("pers",names(totals0[i]))){
      totals0[[i]] <- merge(totals0[[i]][,mget(c(byVars,"Freq"))],
                            data0[,list(FreqPop=sum(get(numberPop))),by=c(byVars)],
                            by=c(byVars))
    }else{
      totals0[[i]] <- merge(totals0[[i]][,mget(c(byVars,"Freq"))],
                            data0[firstPersonInHousehold==TRUE,list(FreqPop=sum(get(numberPop))),by=c(byVars)],
                            by=c(byVars))
    }
  }
  
  return(totals0)
}
# update totals0 with population totals
updateTotals <- function(totals0,data0,hhid,numberPop="weight_choose"){
  firstPersonInHousehold <- NULL
  cn <- colnames(data0)
  TF <- grepl("pers",names(totals0))
  totals0[TF] <- lapply(totals0[TF], function(t0){
    byVars <- intersect(colnames(t0),cn)
    data0[,list(FreqPop=sum(get(numberPop))),by=c(byVars)][t0[,mget(c(byVars,"Freq"))],
                          on=c(byVars)]
  })
  totals0[!TF] <- lapply(totals0[!TF], function(t0){
    byVars <- intersect(colnames(t0),cn)
    data0[firstPersonInHousehold==TRUE,list(FreqPop=sum(get(numberPop))),by=c(byVars)][t0[,mget(c(byVars,"Freq"))],on=c(byVars)]
    })
  return(totals0)
}
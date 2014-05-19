syRec <- function(x, basisvar, hid, newvar, weights,
                  size=sum(weights)/nrow(x), strata){
  if(!is.factor(x[,newvar])) {
    stop("newvar is not a factor. synthetic reconstruction works only with factors.")
  }
  pop <- sampHH(x, sizefactor=size, 
                      hid="db030", strata="db040", 
                      hsize="hsize")[,c(basisvar,hid),with=FALSE]
  ## simulate new variable  
  stat <- as.numeric(levels(eusilcS[,newvar]))
  ## get conditional probablities:
  vars <- c(basisvar,newvar)
  TAB1 <- prop.table(table(eusilcS[,vars]))
  ## simulate new variable:
#  pop <- data.frame(age=eusilcS[eusilcS$db030==hh1,"age"],
#                    sex=eusilcS[eusilcS$db030==hh1,"rb090"])
  samp <- function(x){
    if(all(TAB1[as.numeric(x[1]),as.numeric(x[2]),]==0)){
      status <- NA
    } else {
      status <- sample(stat, 1, prob=TAB1[as.numeric(x[1]),as.numeric(x[2]),])
    }
    as.numeric(status)
  }
  pop <- data.frame(pop)
  stat <- as.numeric(levels(eusilcS[,newvar]))
  pop$status <- aaply(pop, 1, .expand=FALSE, .fun=function(x){
    if(all(TAB1[as.numeric(x[1]),as.numeric(x[2]),]==0)){
      status <- NA
    } else {
      status <- sample(stat, 1, 
                       prob=TAB1[as.numeric(x[1]),as.numeric(x[2]),])
    }
    as.numeric(status)
  })
  pop   
}

install_bitbucket("simPopulation2", username="bernhard_da",  auth_user
                  = "matthias-da", password = "Hdjexly1")

library(synthPop)
data(eusilcS)
selectHHs <- sampHH(eusilcS, sizefactor=sum(weights)/nrow(eusilcS), 
                         hid="db030", strata="db040", hsize="hsize")

syRec(eusilcS, basisvar=c("age","rb090"), hid="db030",
      newvar="pl030", weights="rb050", size=0.01, strata="db040")



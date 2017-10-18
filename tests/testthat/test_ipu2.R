library(simPop)
data(eusilcS)
setDT(eusilcS)
eusilcS <- eusilcS[, list(db030,hsize,db040,age,rb090,netIncome,db090,rb050)]

## some recoding
# generate age groups
eusilcS[age<0, age:=0]
eusilcS[,agegroup:=floor(age/10)]
# some recoding of netIncome for reasons of simplicity
eusilcS[is.na(netIncome), netIncome:=0] 
eusilcS[netIncome<0, netIncome:=0] 
# set hsize to 1,...,5+
eusilcS[hsize>=5, hsize:=5] 

## example for base weights assuming a simple random sample of households stratified per region
eusilcS[, regSamp:=.N, by=db040]
eusilcS[, regPop:=sum(rb050), by=db040]
eusilcS[, baseWeight:=regPop/regSamp]

## constraints on person level
# age 
conP1 <- xtabs(V1 ~ agegroup, data=eusilcS[,sum(rb050),by=agegroup])
# gender by region
conP2 <- xtabs(V1 ~ rb090+db040, data=eusilcS[,sum(rb050),by=list(rb090,db040)])
# personal net income by gender
conP3 <- xtabs(V1 ~ rb090, data=eusilcS[,sum(rb050*netIncome),by=rb090])
## constraints on household level
conH1 <- xtabs(V1 ~ hsize+db040, data=eusilcS[!duplicated(db030),sum(rb050),list(hsize,db040)])
# array of convergence limits for conH1
epsH1 <- conH1
epsH1[as.character(1:4),] <- 0.005
epsH1["5",] <- 0.2

test_that("ipu2 test 1",{
  
  # without array epsP1
  calibweights1 <- ipu2(eusilcS, hid = "db030", 
                        conP = list(conP1,conP2,netIncome=conP3), 
                        conH = list(conH1), 
                        epsP = list(1e-06,1e-06,1e-03),
                        epsH = 0.01,  
                        bound = NULL, verbose=TRUE,  maxIter = 200)
})
test_that("ipu2 test 2",{
  # with array epsP1, base weights and bound
  calibweights2 <- ipu2(eusilcS, hid = "db030", 
                        conP = list(conP1,conP2), 
                        conH = list(conH1), 
                        epsP = 1e-06,
                        epsH = list(epsH1),  
                        w="baseWeight",
                        bound = 4, verbose=TRUE,  maxIter = 200)
})
test_that("ipu2 test 3",{
  # with allPthenH=TRUE
  calibweights2 <- ipu2(eusilcS, hid = "db030", 
                        conP = list(conP1,conP2), 
                        conH = list(conH1), 
                        epsP = 1e-06,
                        epsH = list(epsH1),  
                        w="baseWeight",allPthenH = TRUE,
                        bound = 4, verbose=TRUE,  maxIter = 200)
  
})
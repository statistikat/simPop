################################################
# test calibPop functionality
#
if(Sys.info()["user"]%in%c("kowarik","gussenbauer")){
  

context("calibPop")
library(simPop)

data(eusilcS) # load sample data
data(eusilcP) # population data
#Reduce data to 2 states for computation times
eusilcP <- eusilcP[eusilcP$region%in%c("Vorarlberg","Tyrol"),]
eusilcS <- eusilcS[eusilcS$db040%in%c("Vorarlberg","Tyrol"),]
eusilcP$region <- as.factor(as.character(eusilcP$region))
eusilcS$db040 <- as.factor(as.character(eusilcS$db040))


inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"), method="multinom", nr_cpus=1)

# add margins
margins <- as.data.frame(
  xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$gender + eusilcP$citizenship))
colnames(margins) <- c("db040", "rb090", "pb220a", "freq")
simPop <- addKnownMargins(simPop, margins)

test_that("Test CalibPop - Comparison memory intense vs. memory saving method",{
  simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.01)
  simPop_adj2 <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.01,memory=TRUE)
  expect_equal(simPop_adj2@table[,sum(N)],simPop_adj@table[,sum(N)])
})

test_that("Test CalibPop - check temp.factor",{
  #simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.001,memory=TRUE,choose.temp.factor = .9)
  simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,choose.temp.factor = .5)
  expect_true(abs(simPop_adj@table[,sum(N)]-sum(margins$freq))<1)
})

test_that("Test CalibPop - check sizefactor",{
  simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5)  
  expect_true(abs(simPop_adj@table[,sum(N)]-sum(margins$freq))<1)
})

test_that("Test CalibPop - check scale.redraw",{
  simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5,scale.redraw = .2)
  #simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5,scale.redraw = .8)
  expect_true(abs(simPop_adj@table[,sum(N)]-sum(margins$freq))<1)
})

test_that("Test CalibPop - check observe.break",{
  simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5,observe.break = 0)
  expect_true(abs(simPop_adj@table[,sum(N)]-sum(margins$freq))<1)
  #simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5,observe.break = .3)
})

test_that("Test CalibPop - check observe.times",{
  simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5,observe.times=10)
  #simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5,observe.times=0)
  #simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5,observe.times=10,observe.break = 0.01)
  #simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5,observe.times=10,observe.break = .5)
  expect_true(abs(simPop_adj@table[,sum(N)]-sum(margins$freq))<1)
  simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.1,memory=TRUE,sizefactor = 5,observe.times=10,observe.break = .5)
  expect_true(abs(simPop_adj@table[,sum(N)]-sum(margins$freq))<1)
})

}
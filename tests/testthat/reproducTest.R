#Test for reproducibility
library(simPop)
test_that("Reproduc test",{
data(eusilcS)
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
## in the following, nr_cpus are selected automatically
simPop <- simStructure(data=inp,
                       method="direct", basicHHvars=c("age", "rb090"),seed=10)
simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"),
                         method="multinom", nr_cpus=1,seed=10)
x1 <- simPop@pop@data
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
## in the following, nr_cpus are selected automatically
simPop <- simStructure(data=inp,
                       method="direct", basicHHvars=c("age", "rb090"),seed=10)
simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"),
                         method="multinom", nr_cpus=1,seed=10)
x2 <- simPop@pop@data

if(!identical(x1,x2)){
  stop("setting seed does not work!\n")
}

##Parallel
simPop <- simStructure(data=inp,
                       method="direct", basicHHvars=c("age", "rb090"),seed=10)
simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"),
                         method="multinom", nr_cpus=2,seed=10)
x1 <- simPop@pop@data


simPop <- simStructure(data=inp,
                       method="direct", basicHHvars=c("age", "rb090"),seed=10)
simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"),
                         method="multinom", nr_cpus=2,seed=10)
x2 <- simPop@pop@data
# mclapply can handle parallel seeds, but does not work under Windows.
if(Sys.info()["sysname"] != "Windows"){
  expect_identical(x1,x2,info="setting seed in parallel mode does not work!")
}
})
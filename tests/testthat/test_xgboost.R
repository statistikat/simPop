#Test for xgboost implementation
library(simPop)

test_that("xgboost integration tests",{
  
  data(eusilcS)
  
  inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
  
  ## in the following, nr_cpus are selected automatically
  simPop <- simStructure(data=inp,
                         method="direct", 
                         basicHHvars=c("age", "rb090"),
                         seed=10)
  
  simPop <- simCategorical(simPop,
                           additional=c("pl030", "pb220a"),
                           method="xgboost",
                           nr_cpus=1,
                           seed=10)
  
  eusilcM <- simContinuous(simPop, additional="netIncome",
                           method = "xgboost",
                           regModel = ~rb090+hsize+pl030+pb220a,
                           upper=200000,
                           equidist=FALSE, 
                           nr_cpus=1)
  
  x1 <- eusilcM@pop@data
  
  
  expect_gt(nrow(x1), 0,
            "TODO: failure message")
  
  
  simPop <- simStructure(data=inp,
                         method="direct", basicHHvars=c("age", "rb090"),
                         seed=10)
  
  expect_output({simPop <- simCategorical(simPop,
                                          additional=c("pl030", "pb220a"),
                                          method="xgboost", nr_cpus=1, seed=10, verbose = TRUE)},
                "are running xgboost",
                "TODO: failure message")
  
  expect_output({simPop <- simContinuous(simPop, additional="netIncome",
                                         method = "xgboost",
                                         regModel = ~rb090+hsize+pl030+pb220a,
                                         upper=200000,
                                         equidist=FALSE, 
                                         nr_cpus=1)},
                "are running xgboost",
                "TODO: failure message")
    
})

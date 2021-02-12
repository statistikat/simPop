#Test for xgboost implementation
library(simPop)

# test_that("xgboost integration tests",{
#   
#   data(eusilcS)
#   
#   inp <- specifyInput(data=eusilcS,
#                       hhid="db030",
#                       hhsize="hsize",
#                       strata="db040",
#                       weight="db090")
#   
#   simPop <- simStructure(data=inp,
#                          method="direct", 
#                          basicHHvars=c("age", "rb090"))
#   
#   simPop <- simCategorical(simPop,
#                            additional=c("pl030", "pb220a"),
#                            method="xgboost",
#                            nr_cpus=1)
#   
#   eusilcM <- simContinuous(simPop, additional="netIncome",
#                            method = "xgboost",
#                            regModel = "available",
#                            nr_cpus=1,
#                            by = "none")
#   
#   x1 <- eusilcM@pop@data
#   
#   
#   expect_gt(nrow(x1), 0,
#             "TODO: failure message")
#   
#   
#   simPop <- simStructure(data=inp,
#                          method="direct", basicHHvars=c("age", "rb090"),
#                          seed=10)
#   
#   expect_output({simPop <- simCategorical(simPop,
#                                           additional=c("pl030", "pb220a"),
#                                           method="xgboost",
#                                           nr_cpus=1,
#                                           verbose = TRUE)},
#                 "are running xgboost",
#                 "TODO: failure message")
#   
#   expect_output({simPop <- simContinuous(simPop, additional="netIncome",
#                                          method = "xgboost",
#                                          regModel = "available",
#                                          verbose = TRUE,
#                                          nr_cpus = 1,
#                                          by = "none")},
#                 "are running xgboost",
#                 "TODO: failure message")
#   
#   # TODO: test, regModel, optional_params
# 
# })


test_that("xgboost integration tests",{
  
  data("eusilc13puf")
  data("totalsRGtab")
  
  eusilc13puf$age <- as.numeric(eusilc13puf$age)
  colnames(totalsRGtab) <- c("AT11", "AT21", "AT12", "AT32", "AT22", "AT33", "AT31", "AT13", "AT34")
  totalsRGtab <- totalsRGtab / 100
  
  
  
  inp <- specifyInput(data=eusilc13puf,
                      hhid="db030",
                      hhsize="hsize",
                      strata="db040",
                      weight="rb050")
  
  addWeights(inp) <- calibSample(inp, totalsRGtab) 
  
  simPop <- simStructure(data=inp,
                         method="direct", 
                         basicHHvars=c("age", "rb090"))
  
  simPop <- simCategorical(simPop,
                           additional=c("pl031", "pb220a"),
                           method="xgboost",
                           nr_cpus=1)
  
  eusilcM <- simContinuous(simPop, 
                           additional="hgrossminus",
                           method = "xgboost",
                           regModel = "available",
                           # log = T,
                           nr_cpus=1,
                           by = "db040")
  
  x1 <- eusilcM@pop@data
  
  
  expect_gt(nrow(x1), 0,
            "TODO: failure message")
  
  
  simPop <- simStructure(data=inp,
                         method="direct",
                         basicHHvars=c("age", "rb090"),
                         seed=10)

  expect_output({simPop <- simCategorical(simPop,
                                          additional=c("pl031", "pb220a"),
                                          method="xgboost",
                                          nr_cpus=1,
                                          verbose = TRUE)},
                "are running xgboost",
                "TODO: failure message")

  # TODO: test with by = "none"
  expect_output({simPop <- simContinuous(simPop,
                                         additional="hgrossminus",
                                         method = "xgboost",
                                         regModel = "available",
                                         verbose = TRUE,
                                         nr_cpus = 1,
                                         by = "db040")},
                "are running xgboost",
                "TODO: failure message")
  
  # TODO: test, regModel, optional_params
  
})

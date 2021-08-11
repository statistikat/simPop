#Test for xgboost implementation
library(simPop)

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
                           nr_cpus = 1)
  
  simPop <- simContinuous(simPop, 
                           additional="hgrossminus",
                           method = "xgboost",
                           regModel = "available",
                           by = "db040",
                           log = FALSE,
                           alpha = NULL,
                           nr_cpus = 1)

  expect_gt(nrow(simPop@pop@data), 0,
            "Expected generated synthetic population to have some rows")
  
  
  simPop <- simStructure(data=inp,
                         method="direct",
                         basicHHvars=c("age", "rb090"),
                         seed=10)

  expect_output({simPop <- simCategorical(simPop,
                                          additional=c("pl031", "pb220a"),
                                          method="xgboost",
                                          verbose = TRUE,
                                          nr_cpus = 1)},
                "are running xgboost",
                "Expected commandline output when verbose is true")

  expect_output({simPop <- simContinuous(simPop,
                                         additional="hgrossminus",
                                         method = "xgboost",
                                         regModel = "available",
                                         verbose = TRUE,
                                         by = "db040",
                                         log = FALSE,
                                         alpha = NULL,
                                         nr_cpus = 1)},
                "are running xgboost",
                "Expected commandline output when verbose is true")
})

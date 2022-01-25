#Test for xgboost implementation
library(simPop)

# xgboost integration tests",{
  
  data(eusilcS) # load sample data
  
  ## approx. 20 seconds computation time
  inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
  ## in the following, nr_cpus are selected automatically
  simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
  
  grid <- expand.grid(nrounds = c(10, 50),
                      max_depth = 32,
                      eta = c(0.2, 0.3, 0.5),
                      eval_metric = "mlogloss",
                      stringsAsFactors = F)
  
  simPop <- crossValidation(simPop, additionals=c("pl030", "pb220a"),
                            nr_cpus=1,
                            verbose = TRUE, hyper_param_grid = grid)
  
  expect_true(nrow(simPop@pop@data)> 0,
            "Expected generated synthetic population to have some rows")
#
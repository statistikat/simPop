#Test for reproducibility
library(simPop)

test_that("Keras test", {
  data(eusilcS)
  
  inp <-
    specifyInput(
      data = eusilcS,
      hhid = "db030",
      hhsize = "hsize",
      strata = "db040",
      weight = "db090"
    )
  
  ## in the following, nr_cpus are selected automatically
  simPop <- simStructure(
    data = inp,
    method = "direct",
    basicHHvars = c("age", "rb090"),
    seed = 10
  )
  simPop <- simCategorical(
    simPop,
    additional = c("pl030", "pb220a"),
    method = "keras",
    nr_cpus = 1,
    seed = 10
  )
  
  x1 <- simPop@pop@data
  
  # TODO: test
})
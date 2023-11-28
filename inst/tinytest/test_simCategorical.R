run_tests <- length(strsplit(as.character(packageVersion("simPop")), "[.]")[[1]]) > 3
if(run_tests){
  #Issue 18
  library(simPop)
  data(eusilcS)
  inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", population = TRUE)
  simPop0 <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
  simPop <- simCategorical(simPop0, additional=c("pl030", "pb220a"),
                           method="multinom", nr_cpus=1)
  expect_true(all(c("pl030", "pb220a")%in%colnames(simPop@sample@data)))
  
  simPop <- simCategorical(simPop0, additional=c("pl030", "pb220a"),
                           method="xgboost", nr_cpus=1)
  expect_true(all(c("pl030", "pb220a")%in%colnames(simPop@sample@data)))
  
  simPop <- simCategorical(simPop0, additional=c("pl030", "pb220a"),
                           method="ranger", nr_cpus=1)
  expect_true(all(c("pl030", "pb220a")%in%colnames(simPop@sample@data)))
  
  simPop <- simCategorical(simPop0, additional=c("pl030", "pb220a"),
                           method="distribution", nr_cpus=1)
  expect_true(all(c("pl030", "pb220a")%in%colnames(simPop@sample@data)))
  
  ## with weights
  inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight = "rb050")
  
  simPop0 <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
  simPop <- simCategorical(simPop0, additional=c("pl030", "pb220a"),
                           method="multinom", nr_cpus=1)
  expect_true(all(c("pl030", "pb220a")%in%colnames(simPop@sample@data)))
  
  simPop <- simCategorical(simPop0, additional=c("pl030", "pb220a"),
                           method="xgboost", nr_cpus=1)
  expect_true(all(c("pl030", "pb220a")%in%colnames(simPop@sample@data)))
  
  simPop <- simCategorical(simPop0, additional=c("pl030", "pb220a"),
                           method="ranger", nr_cpus=1)
  expect_true(all(c("pl030", "pb220a")%in%colnames(simPop@sample@data)))
  
  simPop <- simCategorical(simPop0, additional=c("pl030", "pb220a"),
                           method="distribution", nr_cpus=1)
  expect_true(all(c("pl030", "pb220a")%in%colnames(simPop@sample@data)))
}
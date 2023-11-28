run_tests <- length(strsplit(as.character(packageVersion("simPop")), "[.]")[[1]]) > 3
if(run_tests){
  #Test for cforest implementation
  library(simPop)
  
  # cforest integration tests",{
  
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
                           method="cforest",
                           nr_cpus = 1)

  expect_true(nrow(simPop@pop@data)>0,
              "Expected generated synthetic population to have some rows")
  #
}
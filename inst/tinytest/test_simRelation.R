# Tests for simRelation
message("simRelation")
library(simPop)
data(ghanaS) # load sample data
ghanaS[ghanaS$region=="western","nation"] <- ghanaS[ghanaS$region=="western","nation"][1]
samp <- specifyInput(
  data = ghanaS,
  hhid = "hhid",
  strata = "region",
  weight = "weight"
)
ghanaP <- simStructure(
  data = samp,
  method = "direct",
  basicHHvars = c("age", "sex", "relate")
)
class(ghanaP)


# simRelation",{
## long computation time ...
ghanaP <- simRelation(
  simPopObj = ghanaP,
  relation = "relate",
  head = "head",
  additional = c("nation", "ethnic", "religion"), nr_cpus = 1
)
expect_true("nation"%in%colnames(ghanaP@pop@data))
#
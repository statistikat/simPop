context("ipu2")
library(simPop)
data(eusilcS)
setDT(eusilcS)
eusilcS <- eusilcS[, list(db030,hsize,db040,age,rb090,netIncome,db090,rb050)]

## rename columns
setnames(eusilcS, "rb090", "gender")
setnames(eusilcS, "db040", "state")
setnames(eusilcS, "db030", "household")
setnames(eusilcS, "rb050", "weight")

## some recoding
# generate age groups
eusilcS[, agegroup := cut(age, c(-Inf, 10*1:9, Inf), right = FALSE)]
# some recoding of netIncome for reasons of simplicity
eusilcS[is.na(netIncome), netIncome := 0] 
eusilcS[netIncome < 0, netIncome := 0] 
# set hsize to 1,...,5+
eusilcS[, hsize := cut(hsize, c(0:4, Inf), labels = c(1:4, "5+"))]
# treat households as a factor variable
eusilcS[, household := as.factor(household)]

## example for base weights assuming a simple random sample of households stratified per region
eusilcS[, regSamp := .N, by = state]
eusilcS[, regPop := sum(weight), by = state]
eusilcS[, baseWeight := regPop/regSamp]

## constraints on person level
# age 
conP1 <- xtabs(weight ~ agegroup, data = eusilcS)
# gender by region
conP2 <- xtabs(weight ~ gender + state, data = eusilcS)
# personal net income by gender
conP3 <- xtabs(weight*netIncome ~ gender, data = eusilcS)

## constraints on household level
conH1 <- xtabs(weight ~ hsize + state, data = eusilcS, subset = !duplicated(household))

# array of convergence limits for conH1
epsH1 <- conH1
epsH1[1:4,] <- 0.005
epsH1["5+",] <- 0.2





test_that("ipu2 with a numerical variable works as expected", {
  # without array epsP1
  calibweights1 <- ipu2(eusilcS, hid = "household", 
                        conP = list(conP1, conP2, netIncome = conP3), 
                        conH = list(conH1), 
                        epsP = list(1e-06, 1e-06, 1e-03),
                        epsH = 0.01,  
                        bound = NULL, verbose = FALSE,  maxIter = 200)
  expect_true(abs(calibweights1[,sum(calibWeight*netIncome)]-sum(conP3))/sum(conP3)<.01)
  expect_true(all(abs(calibweights1[,sum(calibWeight*netIncome),by=gender]$V1-conP3)/conP3<.01))
  
  err <- max(c(max(abs(xtabs(calibWeight ~ agegroup, data = calibweights1)-conP1)/conP1),
               max(abs(xtabs(calibWeight ~ gender + state, data = calibweights1)-conP2)/conP2),
               max(abs(xtabs(calibWeight ~ hsize + state, data = calibweights1, subset = !duplicated(household))-conH1)/conH1)))
  expect_true(err<.01)
})

test_that("ipu2 works as expected", {
  # with array epsP1, base weights and bound
  calibweights2 <- ipu2(eusilcS, hid = "household", 
                        conP = list(conP1, conP2), 
                        conH = list(conH1), 
                        epsP = 1e-06,
                        epsH = list(epsH1),  
                        w = "baseWeight",
                        bound = 4, verbose = FALSE, maxIter = 200)
  err <- max(c(max(abs(xtabs(calibWeight ~ agegroup, data = calibweights2)-conP1)/conP1),
               max(abs(xtabs(calibWeight ~ gender + state, data = calibweights2)-conP2)/conP2),
               max(abs(xtabs(calibWeight ~ hsize + state, data = calibweights2, subset = !duplicated(household))-conH1)/conH1)))
  expect_true(err<.01)
})

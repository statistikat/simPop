library(parallel)
library(simPopulation)
library(rpart)
library(klaR) # naiveBayes
library(LiblineaR)
library(stringr)
library(microbenchmark)

seed <- 1234
data(eusilcS)   # load sample data

# 0.1s
system.time({
	eusilcP <- simStructure(eusilcS)
})

sapply(list.files("R", full.names=TRUE), source)
source("R/simCategorical-neu.R")

# non-parallel ~ 24secs
system.time({
	eusilcP2 <- simCategorical2(eusilcS, eusilcP, parallel=FALSE)
})
# multinomial regression from package multinom
system.time({
	eusilcP2_1 <- simCategorical2(eusilcS, eusilcP, method="multinom", parallel=TRUE)
})
# recursive partitioning from package rpart
system.time({
	eusilcP2_2 <- simCategorical2(eusilcS, eusilcP, method="rpart", parallel=TRUE)
})
# naivebayes from package klaR
system.time({
	eusilcP2_3 <- simCategorical2(eusilcS, eusilcP, method="naivebayes", parallel=TRUE)
})
# liblinear
system.time({
	eusilcP2_4 <- simCategorical2(eusilcS, eusilcP, method="liblinear", parallel=TRUE)
})

# results
tab <- rbind(
	table(eusilcP2_1$pl030),
	table(eusilcP2_2$pl030),
	table(eusilcP2_3$pl030),
	table(eusilcP2_4$pl030)
)
barplot(tab, beside=TRUE, legend.text=c("multinom","rpart","nbayes","liblinear"), args.legend=list(x="topright"))

# benchmarking
nrruns <- 5
microbenchmark(
	simCategorical2(eusilcS, eusilcP, method="multinom", parallel=TRUE),
	simCategorical2(eusilcS, eusilcP, method="rpart", parallel=TRUE),
	simCategorical2(eusilcS, eusilcP, method="naivebayes", parallel=TRUE),
	simCategorical2(eusilcS, eusilcP, method="liblinear", parallel=TRUE),
	times = 5)



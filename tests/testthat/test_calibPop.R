################################################
# test calibPop functionality
#

library(simPop)

data(eusilcS) # load sample data
data(eusilcP) # population data
## approx. 20 seconds computation time
inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"), method="multinom", nr_cpus=1)

# add margins
margins <- as.data.frame(
  xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$gender + eusilcP$citizenship))
colnames(margins) <- c("db040", "rb090", "pb220a", "freq")
simPop <- addKnownMargins(simPop, margins)

t <- Sys.time()
simPop_adj <- calibPop(simPop, split="db040", temp=1, eps.factor=0.025)
Sys.time()-t
t <- Sys.time()
simPop_adj2 <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.01,memory=TRUE)
Sys.time()-t

simPop_adj
simPop_adj2@table[,sum(N)]

# check choose.temp.factor
simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,choose.temp.factor = 1)
simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,choose.temp.factor = .5)

# check sizefactor
simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,sizefactor = 5)

# check scale.redraw
simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,sizefactor = 5,scale.redraw = .2)
simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,sizefactor = 5,scale.redraw = .8)

# check observe.break
simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,sizefactor = 5,observe.break = 0)
simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,sizefactor = 5,observe.break = .3)

# check observe.times
simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,sizefactor = 5,observe.times=10)

simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,sizefactor = 5,observe.times=0)

simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.001,memory=TRUE,sizefactor = 5,observe.times=10,observe.break = 0.01)

simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.00001,memory=TRUE,sizefactor = 5,observe.times=10,observe.break = .5)
simPop_adj <- calibPop(simPop, split="db040", temp=1, nr_cpus = 1,eps.factor=0.00001,memory=TRUE,sizefactor = 5,observe.times=10,observe.break = .5)





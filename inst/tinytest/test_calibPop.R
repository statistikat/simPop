################################################
# test calibPop functionality
#
message("calibPop")
library(simPop)

data(eusilcS) # load sample data
data(eusilcP) # population data
#Reduce data to 2 states for computation times
eusilcP <- eusilcP[eusilcP$region%in%c("Vorarlberg","Tyrol"),]
eusilcS <- eusilcS[eusilcS$db040%in%c("Vorarlberg","Tyrol"),]
eusilcP$region <- as.factor(as.character(eusilcP$region))
eusilcS$db040 <- as.factor(as.character(eusilcS$db040))


inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"), method="multinom", nr_cpus=1)

# add margins
margins <- as.data.frame(
  xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$region + eusilcP$gender + eusilcP$citizenship))
colnames(margins) <- c("db040", "rb090", "pb220a", "Freq")
simPop <- addKnownMargins(simPop, margins)


# Test CalibPop - check temp.factor",{
epsP.factor <- c(0.1, 0.05, 0.025, 0.01)
for(epsP in epsP.factor){
  simPop_adj <- calibPop(simPop, split="db040", temp=1, epsP.factor=epsP,choose.temp.factor = .5)
  check_table <- merge(simPop_adj@pop@data[,.N,by=.(db040, rb090, pb220a)], margins, by=c("db040", "rb090", "pb220a"))
  expect_true(check_table[,sum(abs(N-Freq))/sum(N)<epsP])
}

simPop_adj <- calibPop(simPop, split="db040", temp=1, epsP.factor=0.025,sizefactor = 5)  
check_table <- merge(simPop_adj@pop@data[,.N,by=.(db040, rb090, pb220a)], margins, by=c("db040", "rb090", "pb220a"))
expect_true(check_table[,sum(abs(N-Freq))/sum(N)<0.025])

# Test CalibPop - check scale.redraw",{
simPop_adj <- calibPop(simPop, split="db040", temp=1, epsP.factor=0.025,sizefactor = 5,scale.redraw = .2)
check_table <- merge(simPop_adj@pop@data[,.N,by=.(db040, rb090, pb220a)], margins, by=c("db040", "rb090", "pb220a"))
expect_true(check_table[,sum(abs(N-Freq))/sum(N)<0.025])
# 

# Test CalibPop - check observe.break",{
simPop_adj <- calibPop(simPop, split="db040", temp=1, epsP.factor=0.025,sizefactor = 5,observe.break = 0)
check_table <- merge(simPop_adj@pop@data[,.N,by=.(db040, rb090, pb220a)], margins, by=c("db040", "rb090", "pb220a"))
expect_true(check_table[,sum(abs(N-Freq))/sum(N)<0.025])
# 

# Test CalibPop - check observe.times",{
simPop_adj <- calibPop(simPop, split="db040", temp=1, epsP.factor=0.025,sizefactor = 5,observe.times=10)
check_table <- merge(simPop_adj@pop@data[,.N,by=.(db040, rb090, pb220a)], margins, by=c("db040", "rb090", "pb220a"))
expect_true(check_table[,sum(abs(N-Freq))/sum(N)<0.025])
simPop_adj <- calibPop(simPop, split="db040", temp=1, epsP.factor=0.025,sizefactor = 5,observe.times=10,observe.break = .5)
check_table <- merge(simPop_adj@pop@data[,.N,by=.(db040, rb090, pb220a)], margins, by=c("db040", "rb090", "pb220a"))
expect_true(check_table[,sum(abs(N-Freq))/sum(N)<0.025])
# 

# Test calibPop - input data.table
pop_data <- pop(simPop)
pop_adj <- calibPop(pop_data, hid = "db030", pid = "pid.simPop", split="db040", 
                       persTables = margins,
                       temp=1, epsP.factor=0.01, choose.temp.factor = .5)
check_table <- merge(pop_adj[,.N,by=.(db040, rb090, pb220a)], margins, by=c("db040", "rb090", "pb220a"))
expect_true(check_table[,sum(abs(N-Freq))/sum(N)<epsP])

# Test calibPop - input data.frame
pop_data <- as.data.frame(pop(simPop))
pop_adj <- calibPop(pop_data, hid = "db030", pid = "pid.simPop", split="db040", 
                    persTables = margins,
                    temp=1, epsP.factor=0.01, choose.temp.factor = .5)
margin_adj <- xtabs(rep(1, nrow(pop_adj)) ~ pop_adj$db040 + pop_adj$rb090 + pop_adj$pb220a)
margin_adj <- as.data.frame(margin_adj)
colnames(margin_adj) <- c("db040", "rb090", "pb220a", "N")
check_table <- merge(margin_adj, margins, by=c("db040", "rb090", "pb220a"))
expect_true(sum(abs(check_table$N-check_table$Freq))/sum(check_table$N)<epsP)

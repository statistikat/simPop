# Just a minimal test instead of examples
#
# Author: alex
###############################################################################
context("simInitSpatial")
library(simPop)
test_that(
  "simInit Spatial Test",
  {
    data(eusilcS)
    data(eusilcP)
    
    # no districts are available in the population, so we have to generate those
    # we randomly assign districts within "region" in the eusilc population data
    # each hh has the same district
    simulate_districts <- function(inp) {
      hhid <- "hid"
      region <- "region"
      
      a <- inp[!duplicated(inp[, hhid]), c(hhid, region)]
      spl <- split(a, a[, region])
      regions <- unique(inp[, region])
      
      tmpres <- lapply(1:length(spl), function(x) {
        codes <- paste(x, 1:sample(3:9, 1), sep = "")
        spl[[x]]$district <-
          sample(codes, nrow(spl[[x]]), replace = TRUE)
        spl[[x]]
      })
      tmpres <- do.call("rbind", tmpres)
      tmpres <- tmpres[, -c(2)]
      out <- merge(
        inp,
        tmpres,
        by.x = c(hhid),
        by.y = hhid,
        all.x = TRUE
      )
      invisible(out)
    }
    
    eusilcP <- data.table(simulate_districts(eusilcP))
    # we generate the input table using the broad region (variable 'region')
    # and the districts, we have generated before.
    #Generate table with household counts by district
    tabHH <-
      eusilcP[!duplicated(hid), .(Freq = .N), by = .(db040 = region, district)]
    setkey(tabHH, db040, district)
    #Generate table with person counts by district
    tabP <- eusilcP[, .(Freq = .N), by = .(db040 = region, district)]
    setkey(tabP, db040, district)
    
    # we generate a synthetic population
    setnames(eusilcP, "region", "db040")
    setnames(eusilcP, "hid", "db030")
    inp <-
      specifyInput(
        data = eusilcP,
        hhid = "db030",
        hhsize = "hsize",
        strata = "db040",
        population = TRUE
      )
    simPopObj <-
      simStructure(
        data = inp,
        method = "direct",
        basicHHvars = c("age", "gender")
      )
    # use only HH counts
    simPopObj1 <-
      simInitSpatial(
        simPopObj,
        additional = "district",
        region = "db040",
        tspatialHH = tabHH,
        tspatialP = NULL,
        nr_cpus = 1
      )
    test1 <-
      merge(tabHH, simPopObj1@pop@data[!duplicated(db030), .N, by = district], by =
              "district")
    expect_false(test1[, max(abs(Freq - N) / Freq)] > 0.01, info =
                   "Test should result in perfect distribution of the households")
    
    simPopObj <-
      simStructure(
        data = inp,
        method = "direct",
        basicHHvars = c("age", "gender")
      )
    # use only P counts
    simPopObj2 <-
      simInitSpatial(
        simPopObj,
        additional = "district",
        region = "db040",
        tspatialHH = NULL,
        tspatialP = tabP,
        nr_cpus = 1
      )
    test2 <-
      merge(tabP, simPopObj2@pop@data[, .N, by = district], by = "district")
    test2[abs(Freq - N) / Freq > 0.05, .N]
    test2[, summary(abs(Freq - N) / Freq)]
    
    expect_false(!test2[abs(Freq - N) / Freq > 0.05, .N] < 10, info = "Test should result in good distribution of the persons")
    
    simPopObj <-
      simStructure(
        data = inp,
        method = "direct",
        basicHHvars = c("age", "gender")
      )
    simPopObj <-
      simCategorical(
        simPopObj,
        additional = "citizenship",
        method = "ranger",
        nr_cpus = 1
      )
    # use P and HH counts
    simPopObj3 <-
      simInitSpatial(
        simPopObj,
        additional = "district",
        region = "db040",
        tspatialHH = tabHH,
        tspatialP = tabP,
        nr_cpus = 1
      )
    test3 <- merge(tabHH[, .(db040, district, FreqH = Freq)],
                   merge(tabP, simPopObj3@pop@data[, .(np = .N, nh = sum(!duplicated(db030))), by =
                                                     district], by = "district"),
                   by = "district")
    
    test3[abs(Freq - np) / Freq > 0.05, .N]
    expect_false(!test3[abs(Freq - np) / Freq > 0.05, .N] < 10,
                 info = "Test should result in good distribution of the persons")
  margins <- eusilcP[, .(freq = .N), by = .(district, gender, citizenship)]
  marginsSynth <-
    simPopObj3@pop@data[, .(freqS = .N), by = .(district, gender, citizenship)]
  simPopObj3 <- addKnownMargins(simPopObj3, margins)
  #simPop_adj <- calibPop(simPopObj3, split="district", temp=10, eps.factor=0.1)
  #marginsSynth2 <- simPop_adj@pop@data[,.(freqSynth=.N),by=.(district,gender,citizenship)]
  #out <- merge(merge(margins,marginsSynth,by=c("district","gender","citizenship")),marginsSynth2,by=c("district","gender","citizenship"))
  })

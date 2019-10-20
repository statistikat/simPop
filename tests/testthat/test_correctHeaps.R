context("heap")
library(simPop)
test_that("correctHeaps",{
  ## create some artificial data
  age <- rlnorm(10000, meanlog=2.466869, sdlog=1.652772)
  age <- round(age[age < 93])
  barplot(table(age))
  
  ## artificially introduce age heaping and correct it:
  # heaps every 5 years
  year5 <- seq(0, max(age), 5)
  age5 <- sample(c(age, age[age %in% year5]))
  cc5 <- rep("darkgrey", length(unique(age)))
  cc5[year5+1] <- "yellow"
  cs <- correctHeaps(age5, heaps="5year", method="lnorm")
  cs <- correctHeaps(age5, heaps="5year", method="norm")
  cs <- correctHeaps(age5, heaps="5year", method="unif")
  
  # heaps every 10 years
  year10 <- seq(0, max(age), 10)
  age10 <- sample(c(age, age[age %in% year10]))
  cc10 <- rep("darkgrey", length(unique(age)))
  cc10[year10+1] <- "yellow"
  cs10 <- correctHeaps(age10, heaps="10year", method="lnorm")
  cs10 <- correctHeaps(age10, heaps="10year", method="norm")
  cs10 <- correctHeaps(age10, heaps="10year", method="unif")
  
  # the first 5 observations should be unchanged
  i1 <- sample(1:length(age10),5)
  i2 <- sample(1:length(age5),5)
  cs10f <- correctHeaps(age10, heaps="10year", method="lnorm", fixed=i1)
  cs5f <- correctHeaps(age5, heaps="5year", method="lnorm", fixed=i2)
  
  
  expect_identical(cs10f[i1],age10[i1])
  expect_identical(cs5f[i2],age5[i2])
})

test_that("correctSingleHeap",{
  ## create some artificial data
  age <- rlnorm(10000, meanlog=2.466869, sdlog=1.652772)
  age <- round(age[age < 93])
  
  ## artificially introduce an age heap for a specific year
  ## and correct it
  age23 <- c(age, rep(23, length=sum(age==23)))
  cc23 <- rep("darkgrey", length(unique(age)))
  cc23[24] <- "yellow"
  cs <- correctSingleHeap(age23, heap=23, before=2, after=3, method="lnorm")
  cs2 <- correctSingleHeap(age23, heap=23, before=5, after=5, method="lnorm")
  cs <- correctSingleHeap(age23, heap=23, before=2, after=3, method="norm")
  cs2 <- correctSingleHeap(age23, heap=23, before=5, after=5, method="norm")
  cs <- correctSingleHeap(age23, heap=23, before=2, after=3, method="unif")
  cs2 <- correctSingleHeap(age23, heap=23, before=5, after=5, method="unif")
  
  # the first 5 observations should be unchanged
  i <- sample(1:length(age23),5)
  csf <- correctSingleHeap(age23, heap=23, before=5, after=5, method="lnorm", fixed=i)
  expect_identical(csf[i],age23[i])
  
})

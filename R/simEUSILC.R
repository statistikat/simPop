simEUSILC <- function(dataS, hid = "db030", wh = "db090",
  wp = "rb050", hsize = NULL, strata = "db040",
  pid = NULL, age = "age", gender = "rb090",
  categorizeAge = TRUE, breaksAge = NULL,
  categorical = c("pl030", "pb220a"),
  income = "netIncome", method = c("multinom", "twostep"),
  breaks = NULL, lower = NULL, upper = NULL,
  equidist = TRUE, probs = NULL, gpd = TRUE,
  threshold = NULL, est = "moments", const = NULL,
  alpha = 0.01, residuals = TRUE,
  components = c("py010n", "py050n", "py090n", "py100n", "py110n", "py120n", "py130n", "py140n"),
  conditional = c(getCatName(income), "pl030"),
  keep = TRUE, maxit = 500, MaxNWts = 1500,
  tol = .Machine$double.eps^0.5, seed) {

  ## initializations
  if ( !missing(seed) ) {
    set.seed(seed)
  }
  if ( !is.character(age) || length(age) != 1 ) {
    stop("'age' must be a character string specifying exactly one column of 'dataS'.\n")
  }
  if ( !is.character(gender) || length(gender) != 1 ) {
    stop("'gender' must be a character string specifying exactly one column of 'dataS'.\n")
  }
  if(!is.character(income) || length(income) != 1) {
    stop("'income' must be a character string specifying exactly one column of 'dataS'\n.")
  }

  ## simulate household structure
  structure <- c(age, gender)
  inp <- specifyInput(data=dataS, hhid=hid, weight=wh, hhsize=hsize, strata=strata, pid=pid)
  synthPop <- simStructure(dataS=inp, basicHHvars=structure)

  dataS <- synthPop@sample@data
  dataP <- synthPop@pop@data

  # construct age categories (if requested)
  categorizeAge <- isTRUE(categorizeAge)
  if ( categorizeAge ) {
    ageCat <- getCatName(age)
    # check breakpoints
    if ( is.null(breaksAge) ) {
      r <- c(range(dataS[[age]], na.rm=TRUE))
      s <- seq(15, 80, 5)
      breaksAge <- c(r[1], s[s > r[1] & s < r[2]], r[2])
    } else if ( !is.numeric(breaksAge) || length(breaksAge) < 2 ) {
      stop("'breaksAge' must be a numeric vector with length >= 2")
    }
    # categorize
    dataS[[ageCat]] <- as.character(cut(dataS[[age]], breaks=breaksAge, include.lowest=TRUE))
    dataP[[ageCat]] <- as.character(cut(dataP[[age]], breaks=breaksAge, include.lowest=TRUE))
    # use age class as predictor instead of age
    structure <- c(ageCat, gender)
  } else ageCat <- NULL
  synthPop@pop@data <- dataP
  synthPop@sample@data <- dataS

  ## simulate additional categorical variables
  basic <- c(structure, if(is.null(hsize)) "hsize" else hsize)
  synthPop <- simCategorical(synthPop, additional=categorical, maxit=maxit, MaxNWts=MaxNWts)

  dataS <- synthPop@sample@data
  dataP <- synthPop@pop@data

  ## simulate income
  basic <- union(basic, categorical)
  method <- match.arg(method)
  useMultinom <- method == "multinom"
  zeros <- TRUE
  # compute income categories
  exclude <- getExclude(dataS[, c(wp, strata, basic, income), with=F])
  if ( length(exclude) > 0 ) {
    incomeS <- dataS[[income]][-exclude]
    wpS <- dataS[[wp]][-exclude]
  } else {
    incomeS <- dataS[[income]]
    wpS <- dataS[[wp]]
  }
  missingBreaks <- is.null(breaks)
  if ( missingBreaks ) {
    if ( is.null(upper) && gpd ) {
      upper <- Inf
    }
    breaks <- getBreaks(incomeS, wpS, zeros, lower, upper, equidist, probs)
  }
  incomeCat <- getCatName(income)
  dataS[[incomeCat]] <- getCat(dataS[[income]], breaks, zeros)
  synthPop@sample@data <- dataS
  synthPop@pop@data <- dataP

  if ( useMultinom ) {
    # multinomial model with random draws from resulting categories
    if ( is.null(threshold) && missingBreaks && (!isTRUE(equidist) || !is.null(probs)) ) {
      threshold <- breaks[length(breaks)-2]
    }
    synthPop <- simContinuous(synthPop, additional=income, zeros=zeros, breaks=breaks,
      gpd=gpd, threshold=threshold, est=est, keep=TRUE, maxit=maxit, MaxNWts=MaxNWts)
  } else {
    # two-step model
    synthPop <- simContinuous(synthPop, additional=income, method="lm", zeros=zeros,
      breaks=breaks, log=TRUE, const=const, alpha=alpha,
      residuals=residuals, maxit=maxit, MaxNWts=MaxNWts, tol=tol)
    dataP <- synthPop@pop@data
    dataP[[incomeCat]] <- getCat(dataP[[income]], breaks, zeros)
    synthPop@pop@data <- dataP
  }

  ## simulate income components
  synthPop <- simComponents(synthPop, total=income, components=components, conditional=conditional)

  # round income components and adjust income
  dataP <- synthPop@pop@data
  dataP[, (components) := lapply(.SD, round, digits=2), .SDcols = components]
  dataP[[income]] <- rowSums(dataP[, components, with=F])

  # return population data
  if ( !keep ) {
    dataP <- dataP[, setdiff(names(dataP), c(ageCat, incomeCat)), with=F]
  }
  synthPop@pop@data <- dataP
  return(invisible(synthPop))
}


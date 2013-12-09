# ----------------------------------------
# Authors: Andreas Alfons and Stefan Kraft
#          Vienna University of Technology
# ----------------------------------------

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
        components = c("py010n", "py050n", "py090n", 
          "py100n", "py110n", "py120n", "py130n", "py140n"), 
        conditional = c(getCatName(income), "pl030"), 
        keep = TRUE, maxit = 500, MaxNWts = 1500, 
        tol = .Machine$double.eps^0.5, seed) {
    
    ## initializations
    if(!missing(seed)) set.seed(seed)  # set seed of random number generator
    if(!is.character(age) || length(age) != 1) { 
        stop("'age' must be a character string ", 
            "specifying exactly one column of 'dataS'")
    }
    if(!is.character(gender) || length(gender) != 1) { 
        stop("'gender' must be a character string ", 
            "specifying exactly one column of 'dataS'")
    }
    if(!is.character(income) || length(income) != 1) { 
        stop("'income' must be a character string ", 
            "specifying exactly one column of 'dataS'")
    }
    
    ## simulate household structure
    structure <- c(age, gender)
    dataP <- simStructure(dataS, hid=hid, w=wh, hsize=hsize, 
        strata=strata, pid=pid, additional=structure, keep=keep)
    
    # construct age categories (if requested)
    categorizeAge <- isTRUE(categorizeAge)
    if(categorizeAge) {
        ageCat <- getCatName(age)
        # check breakpoints
        if(is.null(breaksAge)) {
            r <- c(range(dataS[, age], na.rm=TRUE))
            s <- seq(15, 80, 5)
            breaksAge <- c(r[1], s[s > r[1] & s < r[2]], r[2])
        } else if(!is.numeric(breaksAge) || length(breaksAge) < 2) {
            stop("'breaksAge' must be a numeric vector with length >= 2")
        }
        # categorize
        dataS[, ageCat] <- as.character(cut(dataS[, age], 
                breaks=breaksAge, include.lowest=TRUE))
        dataP[, ageCat] <- as.character(cut(dataP[, age], 
                breaks=breaksAge, include.lowest=TRUE))
        # use age class as predictor instead of age
        structure <- c(ageCat, gender)
    } else ageCat <- NULL
    
    ## simulate additional categorical variables
    basic <- c(structure, if(is.null(hsize)) "hsize" else hsize)
    dataP <- simCategorical(dataS, dataP, w=wp, strata=strata, 
        basic=basic, additional=categorical, maxit=maxit, MaxNWts=MaxNWts)
    
    ## simulate income
    basic <- union(basic, categorical)
    method <- match.arg(method)
    useMultinom <- method == "multinom"
    zeros <- TRUE
    # compute income categories
    exclude <- getExclude(dataS[, c(wp, strata, basic, income)])
    if(length(exclude)) {
        incomeS <- dataS[-exclude, income]
        wpS <- dataS[-exclude, wp]
    } else {
        incomeS <- dataS[, income]
        wpS <- dataS[, wp]
    }
    missingBreaks <- is.null(breaks)
    if(missingBreaks) {
        if(is.null(upper) && gpd) upper <- Inf
        breaks <- getBreaks(incomeS, wpS, zeros, lower, upper, equidist, probs)
    }
    incomeCat <- getCatName(income)
    dataS[, incomeCat] <- getCat(dataS[, income], breaks, zeros)
    if(useMultinom) {
        # multinomial model with random draws from resulting categories
        if(is.null(threshold) && missingBreaks &&  
                (!isTRUE(equidist) || !is.null(probs))) {
            threshold <- breaks[length(breaks)-2]
        }
        dataP <- simContinuous(dataS, dataP, w=wp, strata=strata, 
            basic=basic, additional=income, zeros=zeros, breaks=breaks, 
            gpd=gpd, threshold=threshold, est=est, keep=TRUE, maxit=maxit, 
            MaxNWts=MaxNWts)
    } else {
        # two-step model
        dataP <- simContinuous(dataS, dataP, w=wp, strata=strata, 
            basic=basic, additional=income, method="lm", zeros=zeros, 
            breaks=breaks, log=TRUE, const=const, alpha=alpha, 
            residuals=residuals, maxit=maxit, MaxNWts=MaxNWts, tol=tol)
        dataP[, incomeCat] <- getCat(dataP[, income], breaks, zeros)
    }
    
    ## simulate income components
    dataP <- simComponents(dataS, dataP, w=wp, total=income, 
        components=components, conditional=conditional)
    
    # round income components and adjust income
    dataP[, components] <- round(dataP[, components], digits=2)
    dataP[, income] <- rowSums(dataP[, components])
    
    # return population data
    if(keep) dataP else dataP[, setdiff(names(dataP), c(ageCat, incomeCat))]
}

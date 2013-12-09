# ----------------------------------------
# Authors: Andreas Alfons and Stefan Kraft
#          Vienna University of Technology
# ----------------------------------------

# TODO: further generalization (more than one stratification variable)

simStructure <- function(dataS, hid = "db030", w = "db090", 
        hsize = NULL, strata = "db040", pid = NULL, 
        additional = c("age", "rb090"), 
        method = c("direct", "multinom", "distribution"), 
        keep = TRUE, seed) {
    
    ##### initializations
    if(!missing(seed)) set.seed(seed)  # set seed of random number generator
    if(length(strata) != 1) { 
        stop("currently 'strata' must specify exactly one column of 'dataS'")
    }
    varNames <- c(hid=hid, w=w, hsize=hsize, 
        strata=strata, pid=pid, additional)
    method <- match.arg(method)
    haveHsize <- !is.null(hsize)
    # sample data
    if(all(varNames %in% colnames(dataS))) dataS <- dataS[, varNames]
    else stop("undefined variables selected")
    # order sample (needed later on)
    if(is.null(pid)) dataS <- dataS[order(dataS[, hid]), ]
    else dataS <- dataS[order(dataS[, hid], dataS[, pid]), ]
    # extract variables
    hid <- dataS[, hid]
    w <- dataS[, w]
    if(haveHsize) hsize <- dataS[, hsize] 
    else hsize <- as.integer(table(hid)[as.character(hid)]) 
    strata <- dataS[, strata]
    # preparations of variables
    strata <- as.factor(strata)
    
    ##### set up household structure
    
    # generate variables on household level (indicated by H)
    # these variables have the same value for all household members, hence 
    # we take the first occurence (faster than aggregation with 'unique')
    hfirst <- !duplicated(hid)
    dataH <- data.frame(hsize=as.factor(hsize)[hfirst], strata=strata[hfirst])
    wH <- w[hfirst]
    households <- tableWt(dataH, wH)  # expected number of households
    
    if(method != "direct") {
        ### simulation of household size
        ls <- colnames(households)  # strata levels
        NH <- apply(households, 2, sum)  # number of households by stratum
        if(method == "multinom") {
            empty <- which(households == 0)
            # TODO: allow setting some more arguments
            mod <- suppressWarnings(multinom(hsize~strata, weights=wH, 
                    data=dataH, trace=FALSE))
            newdata <- data.frame(strata=ls)
            rownames(newdata) <- ls
            probs <- as.matrix(predict(mod, newdata=newdata, type="probs"))
#            # experimental: boosting factor for small probabilities
#            fraction <- nrow(dataH)/nrow(dataPH)
#            boost <- 1/sqrt(fraction)
#            probs[probs < fraction] <- pmin(fraction, boost*probs[probs < fraction])
            hsizePH <- unlist(lapply(ls, 
                    function(l) spSample(NH[l], probs[l,])))
        } else if(method == "distribution") {
            hsizePH <- unlist(lapply(ls, 
                    function(l) spSample(NH[l], households[, l])))
        }
        dataPH <- data.frame(hsize=as.factor(hsizePH), 
            strata=factor(rep(ls, times=NH), levels=ls, ordered=is.ordered(strata)))
        households <- tableWt(dataPH)  # recompute number of households
    }
    
    ### simulation of household structure
    
    # all combinations of strata and household size (that do occur
    # in the sample) are stored in data frame 'grid'
    grid <- expand.grid(hsize=rownames(households), 
        strata=colnames(households), stringsAsFactors=FALSE)
    if(method != "multinom") {
        grid <- grid[households != 0,]  # remove those that do not occur
    }
    ncomb <- nrow(grid)
    # indices of households by strata and household size, with combinations 
    # that do not occur in the sample left out (drop = TRUE)
    split <- split(1:nrow(dataH), dataH, drop = (method != "multinom"))
    # to be on the safe side, names are reconstructed from 'grid' 
    # and the list 'split' is sorted according to those names
    nam <- apply(as.matrix(grid), 1, paste, collapse=".")
    split <- split[nam]
    if(method == "multinom") {
        # combinations with strata that do not occur in the sample borrow from 
        # indices of all households in the sample with the same household size
        donor <- grid[empty, 1]
        split[empty] <- split(1:nrow(dataH), dataH$hsize)[donor]        
    }
    # for each stratum, draw from original sample
    numbers <- lapply(1:ncomb, function(i) {
            n <- households[grid[i, 1], grid[i, 2]]
            w <- wH[split[[i]]]
            p <- w / sum(w)  # probability weights
            spSample(n, p)
        })
    
    ### generation of the household structure for the population
    
    # the sampled household numbers in list 'numbers' are transformed to 
    # indices of data frame 'dataS' and stored in vector 'indices'
    hidH <- hid[hfirst]
    indices <- lapply(1:ncomb, function(i) {
            pn <- which(hid %in% hidH[split[[i]]])
            pnM  <- matrix(pn, nrow=as.numeric(grid[i, 1]))
            c(pnM[, numbers[[i]]])
        })
    indices <- unlist(indices)
    
    # new household IDs are created for sampling from population
    if(method == "direct") {
        # household IDs are constructed from the cells in the table of 
        # population households
        upper <- cumsum(households[households != 0])
        lower <- c(1, upper[-ncomb] + 1)
        hidNew <- lapply(1:ncomb, function(i) {
                rep(lower[i]:upper[i], each=as.numeric(grid[i,1]))
            })
    } else {
        # household IDs are constructed from the population household data
        hidNew <- split(1:nrow(dataPH), dataPH, drop=(method == "distribution"))
        hidNew <- lapply(1:ncomb, function(i) {
                rep(hidNew[[i]], each=as.numeric(grid[i,1]))
            })
    }
    hidNew <- unlist(hidNew)
    
    ##### return simulated population household structure
    # it is much faster to generate replications for each variable 
    # individually than for the data frame as a whole
    
    # build command
    if(length(additional)) {
        expr <- paste(additional, "=dataS[, \"", 
            additional, "\"][indices]", sep="")
        expr <- paste(",", expr, collapse="")
    } else expr <- ""
    command <- paste("data.frame(", 
        if(keep) paste(varNames["hid"], "Sample=hid[indices], ", sep=""), 
        varNames["hid"], "=hidNew, ", 
        if(haveHsize) varNames["hsize"] else "hsize", "=hsize[indices], ", 
        varNames["strata"], "=", 
        if(method == "direct") "strata[indices]" else "dataPH$strata[hidNew]", 
        expr, ")", sep="")
    # evaluate command and return result
    dataP <- eval(parse(text=command))
    if(method != "direct") {
        # reorder rows according to constructed population household data
        dataP <- dataP[order(dataP[, varNames["hid"]]),]
        rownames(dataP) <- 1:nrow(dataP)
    }
    dataP
}

# ----------------------------------------
# Authors: Stefan Kraft and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

# handling of negative values is not a problem since function prop.table
# can handle negative values
# both total and components may contain negative values

simComponents <- function(dataS, dataP, w = "rb050", total = "netIncome",
	    components = c("py010n", "py050n", "py090n", 
          "py100n", "py110n", "py120n", "py130n", "py140n"),
        conditional = c(getCatName(total), "pl030"), 
        replaceEmpty = c("sequential", "min"), seed) {
    
	##### initializations
	if(!missing(seed)) set.seed(seed)  # set seed of random number generator
    if(length(total) != 1) {
        stop("currently only one continuous variable can be split at a time")
    }
    if(!isTRUE(length(components) >= 2)) {
        stop("at least two components required")
    }
	varNames <- c(w = w, total=total, components, conditional)
	replaceEmpty <- match.arg(replaceEmpty)
    N <- nrow(dataP)
	
	# check data
    if(all(varNames %in% colnames(dataS))) dataS <- dataS[, varNames]
    else stop("undefined variables in the sample data")
    if(!all(c(total, conditional) %in% colnames(dataP))) {
        stop("undefined variables in the population data")
    }
    
    # exclude observations 
	include <- function(x) !(is.na(x) | x == 0)
    dataS <- dataS[include(dataS[, total]),]
    indPop <- (1:N)[include(dataP[, total])]
    dataPop <- dataP[indPop,]
    
    # data frame dataFrac contains the fractions of the components
    dataFrac <- prop.table(as.matrix(dataS[, components]), 1)
    
    # matrix simFrac stores the simulated fractions
    simFrac <- matrix(ifelse(is.na(dataP[, total]), NA, 0), 
        nrow=N, ncol=length(components))
    colnames(simFrac) <- components
    
    if(length(conditional)) {
        ## class tables
    	tabS <- table(dataS[, conditional])
    	tabP <- table(dataPop[, conditional])
    	if(!identical(dimnames(tabS), dimnames(tabP))) {
            stop("incompatible factor levels in conditioning variables")
        }
        
        ## replace empty combinations in sample by nearest neighbours
        indS <- 1:length(tabS)
        empty <- which(tabS == 0)  # empty cells
        if(length(empty)) {
            donors <- indS[-empty]  # donors
            indTabS <- expand.grid(lapply(dim(tabS), function(k) 1:k)) # indices
            indEmpty <- indTabS[empty, , drop=FALSE]  # indices of empty cells
            indDonors <- indTabS[donors, , drop=FALSE]  # indices of donors
            # reorder donors such that last variable varies fastest
            ord <- do.call("order", indDonors[, names(indDonors), drop=FALSE])
            indDonors <- indDonors[ord, , drop=FALSE]
            donors <- donors[ord]
            # compute replacement cells
            if(replaceEmpty == "sequential" && length(conditional) == 1) {
                replaceEmpty <- "min"
                warning("sequential search of replacement cells not ", 
                    "applicable for only one conditioning variable")
            }
            fun <- switch(replaceEmpty, sequential=seqMinDist, min=minDist)
            replace <- apply(indEmpty, 1, fun, indDonors, donors)
            # replace with donor cell
            indS[empty] <- replace
        }
    	
    	# skip empty combinations in population
    	indP <- which(tabP > 0)
    	indS <- indS[indP]
    	
    	# split the indices of population data by conditioning variables
    	indSplitP <- split(indPop, dataPop[, conditional, drop=FALSE])
    	# split the indices of sample data by conditioning variables
    	indSplitS <- split(1:nrow(dataS), dataS[, conditional, drop=FALSE])
    	# split weights by conditioning variables
    	weights <- split(dataS[, w], dataS[, conditional, drop=FALSE])
    	
    	# sample indices
        indices <- lapply(1:length(indP), 
            function(i) {
                indicesS <- indSplitS[[indS[i]]]
                if(length(indicesS) > 1) {
                    sample(indicesS, size=length(indSplitP[[indP[i]]]), 
                        replace=TRUE, prob=weights[[indS[i]]])
                } else rep.int(indicesS, length(indSplitP[[indP[i]]]))
            })
        
    	# replicate the fractions of the components
    	simFrac[unlist(indSplitP),] <- dataFrac[unlist(indices),]
    } else {
        if(nrow(dataS) > 1) indices <- spSample(length(indPop), dataS[, w])
        else indices <- rep.int(1, length(indPop))
        simFrac[indPop,] <- dataFrac[indices,]
    }
	## this line is excruciatingly slow in R 2.10.0 for large populations
    dataP[, components] <- dataP[, total] * simFrac
	dataP
}

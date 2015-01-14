simComponents <- function(simPopObj, total="netIncome",
  components = c("py010n", "py050n", "py090n", "py100n", "py110n", "py120n", "py130n", "py140n"),
  conditional=c(getCatName(total), "pl030"), replaceEmpty=c("sequential", "min"), seed) {

  ##### initializations
  if ( !missing(seed) ) {
    set.seed(seed)
  }

  dataP <- simPopObj@pop@data
  dataS <- simPopObj@sample@data
  w <- simPopObj@pop@strata

  if ( length(total) != 1 ) {
    stop("currently only one continuous variable can be split at a time\n")
  }
  if ( !isTRUE(length(components) >= 2) ) {
    stop("at least two components are required for this procedure!\n")
  }
  varNames <- c(w=w, total=total, components, conditional)
  replaceEmpty <- match.arg(replaceEmpty)
  N <- nrow(dataP)

  # check data
  if ( all(varNames %in% colnames(dataS)) ) {
    dataS <- dataS[, varNames, with=F]
  } else {
    stop("undefined variables in the sample data!\n")
  }
  if ( !all(c(total, conditional) %in% colnames(dataP)) ) {
    stop("undefined variables in the population data!\n")
  }

  # exclude observations
  include <- function(x) !(is.na(x) | x == 0)
  dataS <- dataS[include(dataS[[total]]),]
  indPop <- (1:N)[include(dataP[[total]])]
  dataPop <- dataP[indPop,]

  # data frame dataFrac contains the fractions of the components
  dataFrac <- prop.table(as.matrix(dataS[, components, with=FALSE]), 1)

  # matrix simFrac stores the simulated fractions
  simFrac <- matrix(ifelse(is.na(dataP[[total]]), NA, 0), nrow=N, ncol=length(components))
  colnames(simFrac) <- components

  if ( length(conditional) > 0 ) {
    ## class tables
    tabS <- table(dataS[, conditional, with=FALSE])
    tabP <- table(dataPop[, conditional, with=FALSE])
    if ( !identical(dimnames(tabS), dimnames(tabP)) ) {
      stop("incompatible factor levels in conditioning variables!\n")
    }

    ## replace empty combinations in sample by nearest neighbours
    indS <- 1:length(tabS)
    empty <- which(tabS == 0)  # empty cells
    if ( length(empty) ) {
      donors <- indS[-empty]  # donors
      indTabS <- expand.grid(lapply(dim(tabS), function(k) 1:k)) # indices
      indEmpty <- indTabS[empty, , drop=FALSE]  # indices of empty cells
      indDonors <- indTabS[donors, , drop=FALSE]  # indices of donors

      # reorder donors such that last variable varies fastest
      ord <- do.call("order", indDonors[, names(indDonors), drop=FALSE])
      indDonors <- indDonors[ord, , drop=FALSE]
      donors <- donors[ord]

      # compute replacement cells
      if ( replaceEmpty == "sequential" && length(conditional) == 1 ) {
        replaceEmpty <- "min"
        warning("sequential search of replacement cells not applicable for only one conditioning variable!\n")
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
    indSplitP <- split(indPop, dataPop[, conditional, with=FALSE])
    # split the indices of sample data by conditioning variables
    indSplitS <- split(1:nrow(dataS), dataS[, conditional, with=FALSE])
    # split weights by conditioning variables
    weights <- split(dataS[[w]], dataS[, conditional, with=FALSE])

    # sample indices
    indices <- lapply(1:length(indP), function(i) {
      indicesS <- indSplitS[[indS[i]]]
      if ( length(indicesS) > 1 ) {
        sample(indicesS, size=length(indSplitP[[indP[i]]]), replace=TRUE, prob=weights[[indS[i]]])
      } else {
        rep.int(indicesS, length(indSplitP[[indP[i]]]))
      }
    })

    # replicate the fractions of the components
    simFrac[unlist(indSplitP),] <- dataFrac[unlist(indices),]
  } else {
    if ( nrow(dataS) > 1 ) {
      indices <- spSample(length(indPop), dataS[[w]])
    } else {
      indices <- rep.int(1, length(indPop))
    }
    simFrac[indPop,] <- dataFrac[indices,]
  }

  out <- dataP[[total]] * simFrac
  for ( i in 1:ncol(out) ) {
    dataP[,colnames(out)[i]] <- out[,i]
  }
  simPopObj@pop@data <- dataP
  invisible(simPopObj)
}

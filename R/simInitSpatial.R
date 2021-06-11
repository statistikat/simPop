#' Generation of smaller regions given an existing spatial variable and a
#' table.
#' 
#' This function allows to manipulate an object of class
#' \code{\linkS4class{simPopObj}} in a way that a new variable containing
#' smaller regions within an already existing broader region is generated. The
#' distribution of the smaller region within the broader region is respected.
#' 
#' The distributional information must be contained in an input table that
#' holds combinations of characteristics of the broader region and the smaller
#' regions as well as population counts (which may be available from a census).
#' 
#' @name simInitSpatial
#' @param simPopObj an object of class \code{\linkS4class{simPopObj}}.
#' @param additional a character vector of length one holding the variable name
#' of the variable containing smaller geographical units. This variable name
#' must be available as a column in input argument \code{tspatial}.
#' @param region a character vector of length one holding the variable name of
#' the broader region. This variable must be available in the input
#' \code{tspatial} as well as in the sample and population slots of input
#' \code{simPopObj}.
#' @param tspatialP a  data.frame (or data.table) containing three columns. The broader region
#' (with the variable name being the same as in input \code{region}, the
#' smaller geographical units (with the variable name being the same as in
#' input \code{additional}) and a third column containing a numeric vector
#' holding counts of persons. This argument or tspatialHH has to be provided.
#' 
#' @param tspatialHH a  data.frame (or data.table) containing three columns. The broader region
#' (with the variable name being the same as in input \code{region}, the
#' smaller geographical units (with the variable name being the same as in
#' input \code{additional}) and a third column containing a numeric vector
#' holding counts of households. This argument or tspatialP has to be provided.
#'
#' @param eps relative deviation of person counts if person and household counts are provided
#' @param maxIter maximum number of iteration for adjustment
#' if person and household counts are provided 
#' @param nr_cpus if specified, an integer number defining the number of cpus
#' that should be used for parallel processing.
#' @param seed optional; an integer value to be used as the seed of the random
#' number generator, or an integer vector containing the state of the random
#' number generator to be restored.
#' @param verbose TRUE/FALSE if some information should be shown during the process
#' @return An object of class \code{\linkS4class{simPopObj}} with an additional
#' variable in the synthetic population slot.
#' @author Bernhard Meindl and Alexander Kowarik
#' @references 
#' M. Templ, B. Meindl, A. Kowarik, A. Alfons, O. Dupriez (2017) Simulation of Synthetic Populations for Survey Data Considering Auxiliary
#' Information. \emph{Journal of Statistical Survey}, \strong{79} (10), 1--38. \doi{10.18637/jss.v079.i10}
#' @keywords manip
#' @export
#' @examples
#' data(eusilcS)
#' data(eusilcP)
#' library(data.table)
#' 
#' # no districts are available in the population, so we have to generate those
#' # we randomly assign districts within "region" in the eusilc population data
#' # each hh has the same district
#' simulate_districts <- function(inp) {
#'   hhid <- "hid"
#'   region <- "region"
#' 
#'   a <- inp[!duplicated(inp[,hhid]),c(hhid, region)]
#'   spl <- split(a, a[,region])
#'   regions <- unique(inp[,region])
#' 
#'   tmpres <- lapply(1:length(spl), function(x) {
#'     codes <- paste(x, 1:sample(3:9,1), sep="")
#'     spl[[x]]$district <- sample(codes, nrow(spl[[x]]), replace=TRUE)
#'     spl[[x]]
#'   })
#'   tmpres <- do.call("rbind", tmpres)
#'   tmpres <- tmpres[,-c(2)]
#'   out <- merge(inp, tmpres, by.x=c(hhid), by.y=hhid, all.x=TRUE)
#'   invisible(out)
#' }
#' 
#' eusilcP <- data.table(simulate_districts(eusilcP))
#' # we generate the input table using the broad region (variable 'region')
#' # and the districts, we have generated before.
#' #Generate table with household counts by district
#' tabHH <- eusilcP[!duplicated(hid),.(Freq=.N),by=.(db040=region,district)]
#' setkey(tabHH,db040,district)
#' #Generate table with person counts by district
#' tabP <- eusilcP[,.(Freq=.N),by=.(db040=region,district)]
#' setkey(tabP,db040,district)
#' 
#' # we generate a synthetic population
#' setnames(eusilcP,"region","db040")
#' setnames(eusilcP,"hid","db030")
#' inp <- specifyInput(data=eusilcP, hhid="db030", hhsize="hsize", strata="db040",population=TRUE)
#' simPopObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "gender"))
#' \donttest{
#' # use only HH counts
#' simPopObj1 <- simInitSpatial(simPopObj, additional="district", region="db040", tspatialHH=tabHH,
#' tspatialP=NULL, nr_cpus=1)
#' 
#' # use only P counts
#' simPopObj2 <- simInitSpatial(simPopObj, additional="district", region="db040", tspatialHH=NULL,
#' tspatialP=tabP, nr_cpus = 1)
#' 
#' # use P and HH counts
#' simPopObj3 <- simInitSpatial(simPopObj, additional="district", region="db040", tspatialHH=tabHH,
#' tspatialP=tabP, nr_cpus = 1)
#' }
#' 
simInitSpatial <- function(simPopObj, additional, region, tspatialP=NULL,tspatialHH=NULL, 
    eps=0.05, maxIter=100, nr_cpus= NULL, seed=1, verbose = FALSE) {
  # simplified version of generateValues_distribution
  generateValues_spatial <- function(dataTable, dataPop, params) {
    # if there is only one subregion in the region
    if( nrow(dataTable) == 1) {
      return(rep(dataTable[[params$additional]],nrow(dataPop)))      
    }
    #give a new name to hhid
    setnames(dataPop,params$predNames,"hhidtmp")
    if(!"freqP"%in%colnames(dataTable)){#Case1: only HH counts provided
      subregion <- sample(rep(dataTable[[params$additional]],times=dataTable[,ceiling(freqPopHH*freqH/sum(freqH))])
      )[1:sum(!duplicated(dataPop[[1]]))]
      dataPop[!duplicated(hhidtmp),subregion:=subregion]
      dataPop[,subregion:=head(subregion,1),by=hhidtmp]
      return(dataPop[["subregion"]])
    
    }else if(!"freqH"%in%colnames(dataTable)){#Case2: only P counts provided
      setnames(dataTable,colnames(dataTable)[2],"subregion")
      dataPop[,subregion:=sample(dataTable[["subregion"]],size=1,prob=dataTable[["freqP"]],
              replace=TRUE),by=hhidtmp]
      dataPop[,.N,by=subregion]
      dataPop[,hsize:=.N,by=hhidtmp]
      meanHH <- dataPop[,mean(hsize)]
      curTab <- merge(dataTable[,list(freqP,subregion)],
          dataPop[,list(nP=.N),by=subregion],by="subregion")
      curTab[,diff:=nP-freqP]
      curTab[,diffp:=diff/freqP]
      setkey(curTab,diffp)
      if(any(abs(curTab[["diffp"]]>eps))){
        dataPop[,firstP:=!duplicated(hhidtmp)]
        for(i in 1:maxIter){
          s <- curTab[nrow(curTab),abs(round(diff/meanHH))]
          x1 <- sample(dataPop[subregion==curTab[nrow(curTab),subregion]&firstP,hhidtmp],size=s)
          tmp <- data.table(hhidtmp=x1,subregionNeu=
                  sample(curTab[diff<0,subregion],size=length(x1),prob=curTab[diff<0,-diff],replace=TRUE))
          dataPop <- merge(dataPop,tmp,by="hhidtmp",all.x=TRUE)
          dataPop[!is.na(subregionNeu),subregion:=subregionNeu]
          dataPop[,subregionNeu:=NULL]
          curTab <- merge(dataTable[,list(freqP,subregion)],
              dataPop[,list(nP=.N),by=subregion],by="subregion")
          curTab[,diff:=nP-freqP]
          curTab[,diffp:=diff/freqP]
          setkey(curTab,diffp)
          if(!any(abs(curTab[["diffp"]]>eps))){
            break;#End for loop if we are close enough
          }
        }  
      }
      
      return(dataPop[["subregion"]])
    }else{#Case3: P and HH counts are provided
      #Initialize the subregion based on the HH counts
      setnames(dataTable,colnames(dataTable)[2],"subregion")
      dataPop[,hsize:=.N,by=hhidtmp]
      meanHH <- dataPop[,mean(hsize)]
      subregion <- sample(rep(dataTable[,subregion],times=dataTable[,ceiling(freqPopHH*freqH/sum(freqH))])
      )[1:sum(!duplicated(dataPop[[1]]))]
      dataPop[!duplicated(hhidtmp),subregion:=subregion]
      dataPop[,subregion:=head(subregion,1),by=hhidtmp]
      #current table
      curTab <- merge(dataTable[,list(freqH,freqP,subregion)],
          dataPop[,list(nP=.N),by=subregion],by="subregion")
      curTab[,diff:=nP-freqP]
      curTab[,diffp:=diff/freqP]
      setkey(curTab,diffp)
      
      if(any(abs(curTab[["diffp"]]>eps))){
        dataPop[,firstP:=!duplicated(hhidtmp)]
        for(i in 1:maxIter){
          s <- curTab[1,abs(round(diff/meanHH))]
          x1 <- sample(dataPop[subregion==curTab[1,subregion]&firstP,hhidtmp],size=s,prob=1/dataPop[subregion==curTab[1,subregion]&firstP,hsize],replace=TRUE)
          x2 <- sample(dataPop[subregion==curTab[nrow(curTab),subregion]&firstP,hhidtmp],size=s,prob=dataPop[subregion==curTab[nrow(curTab),subregion]&firstP,hsize],replace=TRUE)
          dataPop[hhidtmp%in%x2,subregion:=curTab[1,subregion]]
          dataPop[hhidtmp%in%x1,subregion:=curTab[nrow(curTab),subregion]]
          curTab <- merge(dataTable[,list(freqH,freqP,subregion)],
              dataPop[,list(nP=.N),by=subregion],by="subregion")
          curTab[,diff:=nP-freqP]
          curTab[,diffp:=diff/freqP]
          setkey(curTab,diffp)
          if(!any(abs(curTab[["diffp"]]>eps))){
            break;#End for loop if we are close enough
          }
        }  
      }
      return(dataPop[,subregion])
    }
  }
  
  x <- NULL
  diffp <- firstP <- freqH <- freqP <- freqPopHH <- NULL
  hhidtmp <- hsize <- nP <- t <- subregionNeu <- NULL
  Ns <- Nt <- fak <- NULL
  dataP <- simPopObj@pop
  dataS <- simPopObj@sample
  data_pop <- dataP@data
  data_sample <- dataS@data
  basic <- simPopObj@basicHHvars

  
  
  if ( length(additional) != 1 ) {
    stop("currently exactly one additional spatial variable can be generated!\n")
  }
  if ( length(region) != 1 ) {
    stop("argument 'region' must be a character vector of length 1!\n")
  }  
  if ( additional %in% colnames(data_pop) ) {
    stop("variable specified in argument 'additional' already exists in the population!\n")
  }
  if (is.null(tspatialHH)&&is.null(tspatialP)){
    stop("The input for simInitSpatial argument tspatial has changed to tspatialP and tspatialHH
   , one of the 2 has to be defined at least!\n")
  }
  if(!is.null(tspatialHH)){
    tspatialHH <- data.table(tspatialHH)
    if(ncol(tspatialHH)!=3) {
      stop("please check input 'tspatialHH'! Each list element must have exactly 3 columns!\n")
    }
    if(!is.numeric(tspatialHH[[3]])){
      stop("the last and third column of input table tspatialHH must contain numeric values!\n")
    }
    
  }
  if(!is.null(tspatialP)){
    tspatialP <- data.table(tspatialP)
    if(!is.numeric(tspatialP[[3]])){
      stop("the last and third column of input table tspatialP must contain numeric values!\n")
    }
    if(ncol(tspatialP)!=3) {
      stop("please check input 'tspatialP'! Each list element must have exactly 3 columns!\n")
    }
  }
  
  # generation of our main table
  if(!is.null(tspatialHH)){
    tab <- tspatialHH
    colnames(tab) <- c(region,additional,"freqH")
    #Adjust input table to the number in the synth. population
    NsynthHH <- data_pop[!duplicated(data_pop[[simPopObj@pop@hhid]]),list(Ns=.N),by=c(region)]
    tab[,Nt:=sum(freqH),by=c(region)]
    tab <- merge(NsynthHH,tab,by=region)
    tab[,fak:=Ns/Nt]
    tab[,freqH:=freqH*fak]
    tab <- tab[,c(region,additional,"freqH"), with = FALSE]
    if(!is.null(tspatialP)){
      tab <- merge(tab,tspatialP,by=c(region,additional),all=TRUE)
      colnames(tab) <- c(region, additional, "freqH","freqP")
      if(any(is.na(tab[["freqH"]]+tab[["freqP"]]))){
        stop("The table with household counts and person counts dont merge without empty cells.")
      }
    }
  }else{
    tab <- tspatialP
    colnames(tab) <- c(region,additional,"freqP")
  }
  if(!is.null(tspatialP)){
    cn <- colnames(tab)
    tab[,Nt:=sum(freqP),by=c(region)]
    NsynthP <- data_pop[,list(Ns=.N),by=c(region)]
    tab <- merge(tab,NsynthP,by=region)
    tab[,fak:=Ns/Nt]
    tab[,freqP:=freqP*fak]
    tab <- tab[,c(cn),with=FALSE]
  }
  if(verbose){
    cat("The table used for generating the new variable has",nrow(tab),"rows:\n")
    print(tab)
  }

  tab <- merge(tab,data_pop[,list(freqPopP=.N,freqPopHH=sum(!duplicated(eval(parse(text=dataP@hhid))))),by=c(region)],by=region,all=TRUE)
  # Check if the input table match to the synthetic population
  if(any(is.na(rowSums(tab[,na.omit(match(c("freqP","freqH","freqPopP","freqPopHH"),colnames(tab))),with=FALSE])))){
    stop("The table with household counts and person counts does not merge with the population\n without empty cells.")
  }
  # TODO: fix the input (relative?!) if the numbers are not exactly the same
  # to the current synthetic population
  


  # list indStrata contains the indices of dataP split by region
  N <- nrow(data_pop)
  indStrata <- split(1:N, data_pop[[region]])
  
  # predictor variables
  predNames <- dataP@hhid  # in spatial case, it can only be the hhid
  
  params <- list()
  params$additional <- additional
  params$predNames <- predNames
  params$eps <- eps
  params$maxIter <- maxIter
  
  pp <- parallelParameters(nr_cpus=nr_cpus, nr_strata=length(levels(data_sample[[region]])))
  parallel <- pp$parallel
  nr_cores <- pp$nr_cores
  have_win <- pp$have_win; rm(pp)
  
  if ( !missing(seed) ) {
    set.seed(seed,"L'Ecuyer")  # set seed of random number generator
  }
  
  if ( parallel ) {
    # windows
    if ( have_win ) {
      cl <- makePSOCKcluster(nr_cores)
      registerDoParallel(cl,cores=nr_cores)
      values <- foreach(x=levels(data_sample[[region]]), .options.snow=list(preschedule=FALSE)) %dopar% {
        generateValues_spatial(
            dataTable=tab[tab[[region]]==x],
            dataPop=data_pop[indStrata[[x]], predNames, with=FALSE], params)
      }
      stopCluster(cl)
    }
    # linux/max
    if ( !have_win ) {
      values <- mclapply(levels(data_sample[[region]]), function(x) {
            generateValues_spatial(
                dataTable=tab[tab[[region]]==x],
                dataPop=data_pop[indStrata[[x]], predNames, with=FALSE], params)
          }, mc.cores=nr_cores)
    }
  } else {
    values <- lapply(levels(data_sample[[region]]), function(x) {
          generateValues_spatial(
              dataTable=tab[tab[[region]]==x],
              dataPop=data_pop[indStrata[[x]], predNames, with=FALSE], params)
        })  
  }
  names(values) <- levels(data_sample[[region]])
  
  ## add new categorical variables to data set and return
  data_pop[,c(additional):=values[unlist(.BY)],by=c(region)]
  
  # check
  simPopObj@pop@data <- data_pop
  invisible(simPopObj)
}

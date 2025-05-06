#' create an object of class 'dataObj' required for further processing
#' 
#' create an standardized input object of class 'dataObj' containing
#' information on weights, household ids, household sizes, person ids and
#' optionally strata. Outputs of this function are typically used in
#' \code{\link{simStructure}}.
#' 
#' @name specifyInput
#' @param data a \code{data.frame} or \code{data.table} featuring sample data.
#' @param hhid character vector of length 1 specifying variable containing
#' household ids within slot \code{data}. If \code{hhid=NULL} a dummy hhid 
#' (\code{hhid.simPop}) will be created to ensure compatability with
#' other methods/functions in this package.
#' @param hhsize character vector of length 1 specifying variable containing
#' household sizes within slot \code{data}. If NULL, household sizes are
#' automatically calculated.
#' @param pid character vector of length 1 specifying variable containing
#' person ids within slot \code{data}. If NULL, person ids are automatically
#' calculated.
#' @param weight character vector of length 1 specifying variable holding
#' sampling weights within slot \code{data}. If NULL dummy weights \code{weights.simPop=1}
#' will be created to ensure compatability with other methods/functions in this package. 
#' @param strata character vector of length 1 specifing variable name within
#' slot \code{data} of variable holding information on strata, e.g. regions or
#' NULL if such variable does not exist.
#' @param population TRUE/FALSE vector of length 1 specifing if the data object is a sample or a population
#' NULL if such variable does not exist.
#' @export
#' @author Bernhard Meindl
#' @references 
#' M. Templ, B. Meindl, A. Kowarik, A. Alfons, O. Dupriez (2017) Simulation of Synthetic Populations for Survey Data Considering Auxiliary
#' Information. \emph{Journal of Statistical Survey}, \strong{79} (10), 1--38. \doi{10.18637/jss.v079.i10}
#' @keywords method
#' @examples
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", weight="rb050", strata="db040")
#' class(inp)
#' inp
specifyInput <- function(data, hhid=NULL, hhsize=NULL, pid=NULL, weight=NULL,
                         strata=NULL, population=FALSE) {
  
  if (missing(hhid)) hhid <- NULL
  if (missing(hhsize)) hhsize <- NULL
  if (missing(pid)) pid <- NULL
  
  hhid.simPop <- weight.simPop <- NULL
  
  if(!(inherits(data,"data.frame") | inherits(data,"data.table"))){
    stop("data must be either a data.frame or data.table")
  }
  data <- as.data.table(data)
  
  ## ---------------------------
  ## check hhid
  if(is.null(hhid)){
    # initialize hhid
    # used by other methods
    hhid <- "hhid.simPop"
    data[,hhid.simPop := .I]
  }
  if ( !inherits(hhid, "character") | length(hhid) != 1 | is.na(match(hhid, colnames(data)))) {
    stop("hhid must be a character defining the variable holding household ids and must be of length 1!\n")
  }
  setkeyv(data, hhid)
  
  ## ---------------------------
  ## check pid
  if ( !is.null(pid) ) {
    if ( !inherits(pid, "character") | length(pid) != 1 | is.na(match(pid, colnames(data)))) {
      stop("pid must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
    
    if(pid==hhid){
      stop("pid and hhid cannot refer to the same variable in data!")
    }
  } else {
    pid <- "pid.simPop"
    data[,c(pid):=1:.N,by=c(hhid)]
  }
  
  ## ---------------------------
  ## check if weight var is defined well for sample data
  if(!population){
    if(is.null(weight)){
      # weight can be unavailable for data not covering the whole population
      weight <- "weight.simPop"
      data[,weight.simPop:=1]
    }
    if ( !inherits(weight, "character") | length(weight) != 1 | is.na(match(weight, colnames(data)))) {
      stop("weight must be a character defining the variable holding sampling weights and must be of length 1!\n")
    }
  }else{ # for population data no weight var is accepted
    if(!is.null(weight))
      stop("weight must not be defined for a population")
  }
  
  # ---------------------------
  # check strata
  if ( !is.null(strata) ) {
    if ( !inherits(strata, "character") | length(strata) != 1 | is.na(match(strata, colnames(data)))) {
      stop("strata must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
    if(!inherits(data[[strata]], "factor")){
      stop(strata," is not a factor variable as needed for a strata variable.")
    }
    
    if(strata %in% c(pid, hhid)){
      # if strata is same as hhid (or pid) create new variable for hhid (or pid)
      # this will otherwise cause an error when using simCategorical() or simContinuous()
      same_vars <- c("pid", "hhid")[c(pid, hhid) == strata]
      message("strata is identical to argument ",same_vars,", is this intended?\nWill make numerical copy of ",same_vars," which will hold the hhid variable")
      new_id <- paste0(same_vars,".simPop")
      set(data, j = new_id, value = as.numeric(data[[strata]]))
      if(hhid==strata){
        hhid <- new_id
      }
      if(pid==strata){
        pid <- new_id
      }
    }
    
  }else{
    # initialize dummy strata used by other methods in package
    strata <- paste(c("DUMMY_STRATA_",sample(c(letters,LETTERS),8,replace=TRUE)), collapse="")
    data[[strata]] <- factor(1)
  }

  # ---------------------------
  # check hhsize
  if ( !is.null(hhsize) ) {
    if ( !inherits(hhsize, "character") | length(hhsize) != 1 | is.na(match(hhsize, colnames(data)))) {
      stop("hhsize must be a character defining the variable holding information on stratas and must be of length 1!\n")
    }
  } else {
    hhsize <- "hhsize.simPop"
    data[,c(hhsize):=.N,by=c(hhid)]
  }
  
  classes <- sapply(data,class)
  for(i in seq_along(classes)){
    if(classes[[i]][1]=="labelled"){ 
      class(data[[i]]) <- classes[[i]][-1]
    }
  }
  invisible(new("dataObj", data=data, hhid=hhid, hhsize=hhsize, pid=pid, weight=weight, strata=strata, ispopulation=population))
}


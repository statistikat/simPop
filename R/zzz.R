.onAttach <- function(lib, pkg) {
  packageStartupMessage(paste("Package simPop",utils::packageVersion("simPop"),"has been loaded!\n"))
  packageStartupMessage("Since simPop does explicit parallelization,\n the number of data.table threads is set to 1.")
  setDTthreads(1)
}
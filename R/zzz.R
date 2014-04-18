.onAttach <- function(lib, pkg) {
  packageStartupMessage(paste("Package synthPop",utils::packageVersion("synthPop"),"has been loaded!\n"))
}
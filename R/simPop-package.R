#' Simulation of Synthetic Populations for Survey Data 
#' Considering Auxiliary Information
#'
#' The production of synthetic datasets has been proposed as a 
#' statistical disclosure control solution to generate public use 
#' files out of protected data, and as a tool to 
#' create ``augmented datasets'' to serve as input for 
#' micro-simulation models. 
#' Synthetic data have become an important instrument 
#' for \emph{ex-ante} assessments of policies' impact. 
#' The performance and acceptability of such a tool relies 
#' heavily on the quality of the synthetic populations, 
#' i.e., on the statistical similarity between the synthetic 
#' and the true population of interest. 
#' 
#' Multiple approaches and tools have been developed to 
#' generate synthetic data. These approaches can be 
#' categorized into three main groups: 
#' synthetic reconstruction, combinatorial optimization, 
#' and model-based generation. 
#' 
#' The package:
#' \pkg{simPop} is a user-friendly R-package based on a modular object-oriented concept. 
#' It provides a highly optimized S4 class implementation 
#' of various methods, including calibration by iterative 
#' proportional fitting and simulated annealing, and 
#' modeling or data fusion by logistic regression. 
#' 
#' The following applications further shows the methods and package:   
#' We firstly demonstrated the use of \pkg{simPop} by creating 
#' a synthetic population of Austria based on the 
#' European Statistics of Income and Living Conditions (Alfons et al., 2011) 
#' including the evaluation of the quality of the generated population. 
#' In this contribution, the mathematical details of functions \code{simStructure}, \code{simCategorical},
#' \code{simContinuous} and \code{simComponents} are given in detail.
#' The disclosure risk of this synthetic population has been evaluated in (Templ and Alfons, 2012) using large-scale simulation studies. 
#' 
#' Employer-employee data were created in Templ and Filzmoser (2014) whereby 
#' the structure of companies and employees are considered.
#' 
#' Finally, the R package \pkg{simPop} is presented in full detail 
#' in Templ et al. (2017). In this paper - the main reference to this work - 
#' all functions and the S4 class 
#' structure of the package are described in detail. For beginners, this paper might be 
#' the starting point to learn about the methods and package.
#'
#'
#' \tabular{ll}{ Package: \tab simPop\cr Type: \tab Package\cr Version: \tab
#' 1.0.0\cr Date: \tab 20017-08-07\cr License: \tab GPL (>= 2) \cr }
#'
#' @name simPop-package
#' @aliases simPop-package simPop
#' @docType package
#' @author Bernhard Meindl, Matthias Templ, Andreas Alfons, Alexander Kowarik,
#' 
#' Maintainer: Matthias Templ <matthias.templ@@gmail.com>
#' @references 
#' M. Templ, B. Meindl, A. Kowarik, A. Alfons, O. Dupriez (2017) Simulation of Synthetic Populations for Survey Data Considering Auxiliary
#' Information. \emph{Journal of Statistical Survey}, \strong{79} (10), 1--38. \doi{10.18637/jss.v079.i10}
#' 
#' A. Alfons, M. Templ (2011) Simulation of close-to-reality population data for household surveys with application to EU-SILC. \emph{Statistical Methods & Applications}, \strong{20} (3), 383--407. doi: 10.1007/s10260-011-0163-2
#' 
#' M. Templ, P. Filzmoser (2014) Simulation and quality of a synthetic close-to-reality employer-employee population. 
#' Journal of Applied Statistics, \strong{41} (5), 1053--1072. \doi{10.1080/02664763.2013.859237} 
#' 
#' M. Templ, A. Alfons (2012) Disclosure Risk of Synthetic Population Data 
#' with Application in the Case of EU-SILC. In J Domingo-Ferrer, E Magkos (eds.), 
#' \emph{Privacy in Statistical Databases}, \strong{6344} of Lecture Notes in Computer Science, 174--186. Springer Verlag, Heidelberg. \doi{10.1007/978-3-642-15838-4_16}
#' 
#' @keywords package
#' @examples
#'
#' ## we use synthetic eusilcS survey sample data 
#' ## included in the package to simulate a population
#' 
#' ## create the structure
#' data(eusilcS)
#' \donttest{
#' ## approx. 20 seconds computation time
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' ## in the following, nr_cpus are selected automatically
#' simPop <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' simPop <- simCategorical(simPop, additional=c("pl030", "pb220a"), method="multinom", nr_cpus=1)
#' simPop
#' class(simPop)
#' regModel = ~rb090+hsize+pl030+pb220a
#'
#' ## multinomial model with random draws
#' eusilcM <- simContinuous(simPop, additional="netIncome",
#'               regModel = regModel,
#'               upper=200000, equidist=FALSE, nr_cpus=1)
#' class(eusilcM)
#' }
#' 
#' ## this is already a basic synthetic population, but
#' ## many other functions in the package might now 
#' ## be used for fine-tuning, adding further variables, 
#' ## evaluating the quality, adding finer geographical details, 
#' ## using different methods, calibrating surveys or populations, etc. 
#' ## -- see Templ et al. (2017) for more details.
#' 
NULL

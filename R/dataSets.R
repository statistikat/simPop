#' Population totals Region times Gender for Austria 2006
#' 
#' Population characteristics Region times Gender from Austria.
#' 
#' 
#' @name totalsRG
#' @aliases totalsRG totalsRGtab
#' @docType data
#' @format totalsRG: A data frame with 18 observations on the following 3
#' variables.  \describe{ \item{list("rb090")}{gender; a factor with levels
#' \code{female} \code{male}} \item{list("db040")}{region; a factor with levels
#' \code{Burgenland} \code{Carinthia} \code{Lower Austria,} \code{Salzburg}
#' \code{Styria} \code{Tyrol} \code{Upper Austria} \code{Vienna}
#' \code{Vorarlberg}} \item{list("Freq")}{totals; a numeric vector} }
#' totalsRGtab: a two-dimensional table holding the same information
#' @source StatCube - statistical data base,
#' \url{http://www.statistik.at/web_de/services/datenbank_superstar/}
#' @keywords datasets
#' @examples
#' 
#' data(totalsRG)
#' totalsRG
#' data(totalsRGtab)
#' totalsRGtab
#' 
NULL


#' Synthetic EU-SILC data
#' 
#' This data set is synthetically generated from real Austrian EU-SILC
#' (European Union Statistics on Income and Living Conditions) data.
#' 
#' The data set is used as population data in some of the examples in package
#' \code{simFrame}.  Note that it is included for illustrative purposes only.
#' It consists of 25 000 households, hence it does not represent the true
#' population sizes of Austria and its regions.
#' 
#' Only a few of the large number of variables in the original survey are
#' included in this example data set.  Some variable names are different from
#' the standardized names used by the statistical agencies, as the latter are
#' rather cryptic codes.  Furthermore, the variables \code{hsize},
#' \code{eqsize}, \code{eqIncome} and \code{age} are not included in the
#' standardized format of EU-SILC data, but have been derived from other
#' variables for convenience.  Moreover, some very sparse income components
#' were not included in the the generation of this synthetic data set. Thus the
#' equivalized household income is computed from the available income
#' components.
#' 
#' @name eusilcP
#' @aliases eusilcP 
#' @docType data
#' @format A \code{data.frame} with 58 654 observations on the following 28
#' variables: \describe{ \item{hid}{integer; the household ID.}
#' \item{region}{factor; the federal state in which the household is
#' located (levels \code{Burgenland}, \code{Carinthia}, \code{Lower Austria},
#' \code{Salzburg}, \code{Styria}, \code{Tyrol}, \code{Upper Austria},
#' \code{Vienna} and \code{Vorarlberg}).} \item{hsize}{integer; the
#' number of persons in the household.} \item{eqsize}{numeric; the
#' equivalized household size according to the modified OECD scale.}
#' \item{eqIncome}{numeric; a simplified version of the equivalized
#' household income.} \item{pid}{integer; the personal ID.}
#' \item{id}{the household ID combined with the personal ID.  The first five
#' digits represent the household ID, the last two digits the personal ID (both
#' with leading zeros).} \item{age}{integer; the person's age.}
#' \item{gender}{factor; the person's gender (levels \code{male} and
#' \code{female}).} \item{ecoStat}{factor; the person's economic status
#' (levels \code{1} = working full time, \code{2} = working part time, \code{3}
#' = unemployed, \code{4} = pupil, student, further training or unpaid work
#' experience or in compulsory military or community service, \code{5} = in
#' retirement or early retirement or has given up business, \code{6} =
#' permanently disabled or/and unfit to work or other inactive person, \code{7}
#' = fulfilling domestic tasks and care responsibilities).}
#' \item{citizenship}{factor; the person's citizenship (levels
#' \code{AT}, \code{EU} and \code{Other}).} \item{py010n}{numeric;
#' employee cash or near cash income (net).} \item{py050n}{numeric;
#' cash benefits or losses from self-employment (net).}
#' \item{py090n}{numeric; unemployment benefits (net).}
#' \item{py100n}{numeric; old-age benefits (net).}
#' \item{py110n}{numeric; survivor's benefits (net).}
#' \item{py120n}{numeric; sickness benefits (net).}
#' \item{py130n}{numeric; disability benefits (net).}
#' \item{py140n}{numeric; education-related allowances (net).}
#' \item{hy040n}{numeric; income from rental of a property or land
#' (net).} \item{hy050n}{numeric; family/children related allowances
#' (net).} \item{hy070n}{numeric; housing allowances (net).}
#' \item{hy080n}{numeric; regular inter-household cash transfer
#' received (net).} \item{hy090n}{numeric; interest, dividends, profit
#' from capital investments in unincorporated business (net).}
#' \item{hy110n}{numeric; income received by people aged under 16
#' (net).} \item{hy130n}{numeric; regular inter-household cash transfer
#' paid (net).} \item{hy145n}{numeric; repayments/receipts for tax
#' adjustment (net).} \item{main}{logical; indicates the main income
#' holder (i.e., the person with the highest income) of each household.} }
#' @references Eurostat (2004) Description of target variables: Cross-sectional
#' and longitudinal. \emph{EU-SILC 065/04}, Eurostat.
#' @source This is a synthetic data set based on Austrian EU-SILC data from
#' 2006.  The original sample was provided by Statistics Austria.
#' @keywords datasets
#' @examples
#' 
#' data(eusilcP)
#' summary(eusilcP)
#' 
NULL





#' Synthetic EU-SILC survey data
#' 
#' This data set is synthetically generated from real Austrian EU-SILC
#' (European Union Statistics on Income and Living Conditions) data.
#' 
#' The data set consists of 4641 households and is used as sample data in some
#' of the examples in package \code{simPopulation}.  Note that it is included
#' for illustrative purposes only.  The sample weights do not reflect the true
#' population sizes of Austria and its regions.  The resulting population data
#' is about 100 times smaller than the real population size to save computation
#' time.
#' 
#' Only a few of the large number of variables in the original survey are
#' included in this example data set.  The variable names are rather cryptic
#' codes, but these are the standardized names used by the statistical
#' agencies.  Furthermore, the variables \code{hsize}, \code{age} and
#' \code{netIncome} are not included in the standardized format of EU-SILC
#' data, but have been derived from other variables for convenience.
#' 
#' @name eusilcS
#' @docType data
#' @format A data frame with 11725 observations on the following 18 variables.
#' \describe{ \item{db030}{integer; the household ID.}
#' \item{hsize}{integer; the number of persons in the household.}
#' \item{db040}{factor; the federal state in which the household is
#' located (levels \code{Burgenland}, \code{Carinthia}, \code{Lower Austria},
#' \code{Salzburg}, \code{Styria}, \code{Tyrol}, \code{Upper Austria},
#' \code{Vienna} and \code{Vorarlberg}).} \item{age}{integer; the
#' person's age.} \item{rb090}{factor; the person's gender (levels
#' \code{male} and \code{female}).} \item{pl030}{factor; the person's
#' economic status (levels \code{1} = working full time, \code{2} = working
#' part time, \code{3} = unemployed, \code{4} = pupil, student, further
#' training or unpaid work experience or in compulsory military or community
#' service, \code{5} = in retirement or early retirement or has given up
#' business, \code{6} = permanently disabled or/and unfit to work or other
#' inactive person, \code{7} = fulfilling domestic tasks and care
#' responsibilities).} \item{pb220a}{factor; the person's citizenship
#' (levels \code{AT}, \code{EU} and \code{Other}).}
#' \item{netIncome}{numeric; the personal net income.}
#' \item{py010n}{numeric; employee cash or near cash income (net).}
#' \item{py050n}{numeric; cash benefits or losses from self-employment
#' (net).} \item{py090n}{numeric; unemployment benefits (net).}
#' \item{py100n}{numeric; old-age benefits (net).}
#' \item{py110n}{numeric; survivor's benefits (net).}
#' \item{py120n}{numeric; sickness benefits (net).}
#' \item{py130n}{numeric; disability benefits (net).}
#' \item{py140n}{numeric; education-related allowances (net).}
#' \item{db090}{numeric; the household sample weights.}
#' \item{rb050}{numeric; the personal sample weights.} }
#' @references Eurostat (2004) Description of target variables: Cross-sectional
#' and longitudinal. \emph{EU-SILC 065/04}, Eurostat.
#' @source This is a synthetic data set based on Austrian EU-SILC data from
#' 2006.  The original sample was provided by Statistics Austria.
#' @keywords datasets
#' @examples
#' 
#' data(eusilcS)
#' summary(eusilcS)
#' 
NULL


#' Synthetic EU-SILC 2013 survey data
#' 
#' This data set is synthetically generated from real Austrian EU-SILC
#' (European Union Statistics on Income and Living Conditions) data 2013.
#' 
#' The data set consists of 5977 households and is used as sample data in some
#' of the examples in package \code{simPop}.  Note that it is included
#' for illustrative purposes only.  The sample weights do not reflect the true
#' population sizes of Austria and its regions.  
#' 
#' 62 variables of the original survey are
#' simulated for this example data set.  The variable names are rather cryptic
#' codes, but these are the standardized names used by the statistical
#' agencies.  Furthermore, the variables \code{hsize}, \code{age} and
#' \code{netIncome} are not included in the standardized format of EU-SILC
#' data, but have been derived from other variables for convenience.
#' 
#' @name eusilc13puf
#' @docType data
#' @format A data frame with 13513 observations on the following 62 variables.
#' \describe{ 
#' \item{db030}{integer; the household ID.}
#' \item{hsize}{integer; the number of persons in the household.}
#' \item{db040}{factor; the federal state in which the household is
#' located (levels \code{Burgenland}, \code{Carinthia}, \code{Lower Austria},
#' \code{Salzburg}, \code{Styria}, \code{Tyrol}, \code{Upper Austria},
#' \code{Vienna} and \code{Vorarlberg}).} 
#' \item{age}{integer; the
#' person's age.} 
#' \item{rb090}{factor; the person's gender (levels
#' \code{male} and \code{female}).} 
#' \item{pid}{personal ID}
#' \item{weight}{sampling weights}
#' \item{pl031}{factor; the person's
#' economic status (levels \code{1} = working full time, \code{2} = working
#' part time, \code{3} = unemployed, \code{4} = pupil, student, further
#' training or unpaid work experience or in compulsory military or community
#' service, \code{5} = in retirement or early retirement or has given up
#' business, \code{6} = permanently disabled or/and unfit to work or other
#' inactive person, \code{7} = fulfilling domestic tasks and care
#' responsibilities).} 
#' \item{pb220a}{factor; the person's citizenship
#' (levels \code{AT}, \code{EU} and \code{Other}).}
#' \item{...}{see Eurostat's code book for further variables}
#' }
#' @references Eurostat (2013) Description of target variables: Cross-sectional
#' and longitudinal.
#' @source This is a synthetic data set based on Austrian EU-SILC data from
#' 2013.  The original sample was provided by Statistics Austria.
#' @keywords datasets
#' @author Matthias Templ
#' @examples
#' data(eusilc13puf)
#' str(eusilc13puf)
NULL


#' Synthetic GLSS survey data
#' 
#' This data set is synthetically generated from real GLSS (Ghana Living
#' Standards Survey) data.
#' 
#' The data set consists of 8700 households and is used as sample data in some
#' of the examples in package \code{simPopulation}.  Note that it is included
#' for illustrative purposes only.  The sample weights do not reflect the true
#' population sizes of Ghana and its regions.  The resulting population data is
#' about 100 times smaller than the real population size to save computation
#' time.
#' 
#' Only some of the variables in the original survey are included in this
#' example data set.  Furthermore, categories are aggregated for certain
#' variables due to the large number of possible outcomes in the original
#' survey data.
#' 
#' @name ghanaS
#' @docType data
#' @format A data frame with 36970 observations on the following 14 variables.
#' \describe{ \item{hhid}{integer; the household ID.}
#' \item{hsize}{integer; the number of persons in the household.}
#' \item{region}{factor; the region in which the household is located
#' (levels \code{western}, \code{central}, \code{greater accra}, \code{volta},
#' \code{eastern}, \code{ashanti}, \code{brong ahafo}, \code{northern},
#' \code{upper east} and \code{upper west}).} \item{clust}{factor; the
#' enumeration area.} \item{age}{integer; the person's age.}
#' \item{sex}{factor; the person's sex (levels \code{male} and
#' \code{female}).} \item{relate}{factor; the relationship with the
#' household head (levels \code{head}, \code{spouse}, \code{child},
#' \code{grandchild}, \code{parent/parentlaw}, \code{son/daughterlaw},
#' \code{other relative}, \code{adopted child}, \code{househelp} and
#' \code{non_relative}).} \item{nation}{factor; the person's
#' nationality (levels \code{ghanaian birth}, \code{ghanaian naturalise},
#' \code{burkinabe}, \code{malian}, \code{nigerian}, \code{ivorian},
#' \code{togolese}, \code{liberian}, \code{other ecowas}, \code{other africa}
#' and \code{other}).} \item{ethnic}{factor; the person's ethnicity
#' (levels \code{akan}, \code{all other tribes}, \code{ewe}, \code{ga-dangbe},
#' \code{grusi}, \code{guan}, \code{gurma}, \code{mande} and
#' \code{mole-dagbani}).} \item{religion}{factor; the person's religion
#' (levels \code{catholic}, \code{anglican}, \code{presbyterian},
#' \code{methodist}, \code{pentecostal}, \code{spiritualist}, \code{other
#' christian}, \code{moslem}, \code{traditional}, \code{no religion} and
#' \code{other}).} \item{highest_degree}{factor; the person's highest
#' degree of education (levels \code{none}, \code{mlsc}, \code{bece},
#' \code{voc/comm}, \code{teacher trng a}, \code{teacher trng b}, \code{gce 'o'
#' level}, \code{ssce}, \code{gce 'a' level}, \code{tech/prof cert},
#' \code{tech/prof dip}, \code{hnd}, \code{bachelor}, \code{masters},
#' \code{doctorate} and \code{other}).} \item{occupation}{factor; the
#' person's occupation (levels \code{armed forces and other security
#' personnel}, \code{clerks}, \code{craft and related trades workers},
#' \code{elementary occupations}, \code{legislators, senior officials and
#' managers}, \code{none}, \code{plant and machine operators and assemblers},
#' \code{professionals}, \code{service workers and shop and market sales
#' workers}, \code{skilled agricultural and fishery workers}, and
#' \code{technicians and associate professionals}).}
#' \item{income}{numeric; the person's annual income.}
#' \item{weight}{numeric; the sample weights.} }
#' @references Ghana Statistical Service (2008) Ghana Living Standards Survey:
#' Report of the fifth round.
#' @source This is a synthetic data set based on GLSS data from 2006.  The
#' original sample was provided by Ghana Statistical Service.
#' @keywords datasets
#' @examples
#' 
#' data(ghanaS)
#' summary(ghanaS)
#' 
NULL





#' Extract and modify variables from population or sample data stored in an
#' object of class \code{\link{simPopObj-class}}.
#' 
#' Using \code{\link{samp}} \code{\link{samp<-}} it is possible to extract or
#' rather modify variables of the sample data within slot \code{data} in slot
#' \code{sample} of the \code{\link{simPopObj-class}}-object. Using
#' \code{\link{pop}} \code{\link{pop<-}} it is possible to extract or rather
#' modify variables of the synthetic population within in slot \code{data} in
#' slot \code{sample} of the \code{\link{simPopObj-class}}-object.
#' 
#' 
#' @name get_set-methods
#' @aliases pop pop<- samp samp<- samp,simPopObj-method pop,simPopObj-method
#' samp<-,simPopObj-method pop<-,simPopObj-method
#' @docType methods
#' @param obj An object of class \code{\link{simPopObj-class}}
#' @param var variable name or index for the variable in slot 'samp' of object
#' with the slot name to be accessed. If \code{NULL}, the entire dataset
#' (sample or population) is returned.
#' @param value Content replacing whatever the variable in slot \code{var} in
#' \code{obj} currently holds.
#' @return Returns an object of class \code{\link{simPopObj-class}} with the
#' appropriate replacement.
#' @author Bernhard Meindl
#' @seealso \code{\link{simPopObj-class}},\code{\link{pop}},
#' \code{\link{pop<-}}, \code{\link{samp<-}}, \code{\link{manageSimPopObj}}
#' @keywords manip methods
#' @examples
#' 
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040",
#' weight="db090")
#' simPopObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' 
#' ## get/set variables in sample-object of simPopObj
#' head(samp(simPopObj, var="age"))
#' samp(simPopObj, var="newVar") <- 1
#' head(samp(simPopObj, var="newVar"))
#' ## deleting is also possible
#' samp(simPopObj, var="newvar") <- NULL
#' head(samp(simPopObj, var="newvar"))
#' ## extract multiple variables
#' head(samp(simPopObj, var=c("db030","db040")))
#' 
#' ## get/set variables in pop-object of simPopObj
#' head(pop(simPopObj, var="age"))
#' pop(simPopObj, var="newVar") <- 1
#' head(pop(simPopObj, var="newVar"))
#' ## deleting is also possible
#' pop(simPopObj, var="newvar") <- NULL
#' head(pop(simPopObj, var="newvar"))
#' ## extract multiple variables
#' head(pop(simPopObj, var=c("db030","db040")))
#' 
NULL





#' Methods to query/set slots from objects of class
#' \code{\linkS4class{simPopObj}}.
#' 
#' The methods allow to query or replace various information of
#' \code{\linkS4class{simPopObj}}-objects.
#' 
#' 
#' @name simPopMethods
#' @aliases sampleObj sampleObj-methods sampleObj<- sampleObj<--methods
#' sampleData popObj popObj-methods popObj<- popObj<--methods popData tableObj
#' tableObj-methods sampleObj,simPopObj-method
#' sampleObj<-,simPopObj,dataObj-method sampleData,simPopObj-method
#' popObj,simPopObj-method popObj<-,simPopObj,dataObj-method
#' popData,simPopObj-method tableObj,simPopObj-method
#' @docType methods
#' @section Methods: \describe{ Below, the methods to access/set various slots
#' of \code{\linkS4class{simPopObj}}-objects are shown.
#' \item{list("sampleObj")}{ returns slot 'sample' of the input object. }\item{
#' with }{ returns slot 'sample' of the input object.
#' }\item{list("signature(object=\"simPopObj\")")}{ returns slot 'sample' of
#' the input object. } \item{list("sampleObj<-")}{ replaces slot 'sample' of
#' the input object with the object \code{value}. }\item{ with }{ replaces slot
#' 'sample' of the input object with the object \code{value}.
#' }\item{list("signature(object=\"simPopObj\", value=\"dataObj\")")}{ replaces
#' slot 'sample' of the input object with the object \code{value}. }
#' \item{list("sampleData")}{ returns slot 'data' of slot 'sample' of the input
#' object. }\item{ with }{ returns slot 'data' of slot 'sample' of the input
#' object. }\item{list("signature(object=\"simPopObj\")")}{ returns slot 'data'
#' of slot 'sample' of the input object. } \item{list("popObj")}{ returns slot
#' 'pop' of the input object. }\item{ with }{ returns slot 'pop' of the input
#' object. }\item{list("signature(object=\"simPopObj\")")}{ returns slot 'pop'
#' of the input object. } \item{list("popObj<-")}{ replaces slot 'pop' of the
#' input object with the object \code{value}. }\item{ with }{ replaces slot
#' 'pop' of the input object with the object \code{value}.
#' }\item{list("signature(object=\"simPopObj\", value=\"dataObj\")")}{ replaces
#' slot 'pop' of the input object with the object \code{value}. }
#' \item{list("popData")}{ returns slot 'data' of slot 'pop' of the input
#' object. }\item{ with }{ returns slot 'data' of slot 'pop' of the input
#' object. }\item{list("signature(object=\"simPopObj\")")}{ returns slot 'data'
#' of slot 'pop' of the input object. } \item{list("tableObj")}{ returns slot
#' 'table' of the input object. }\item{ with }{ returns slot 'table' of the
#' input object. }\item{list("signature(object=\"simPopObj\")")}{ returns slot
#' 'table' of the input object. } }
#' @keywords methods
#' @examples
#' 
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' inpObj <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' class(inpObj)
#' str(sampleObj(inpObj))
#' 
#' inp2 <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="pb220a", weight="db090")
#' sampleObj(inpObj) <- inp2
#' str(sampleObj(inpObj))
#' 
NULL





#' Class \code{"simPopObj"}
#' 
#' An object that is used throughout the package containing information on the
#' sample (in slot \code{sample}), the population (slot \code{pop}) and
#' optionally some margins in form of a table (slot \code{table}).
#' 
#' 
#' @name simPopObj-class
#' @aliases simPopObj-class show,simPopObj-method
#' @docType class
#' @section Objects from the Class: Objects are automatically created in
#' function \code{\link{simStructure}}.
#' @author Bernhard Meindl and Matthias Templ
#' @seealso \code{\linkS4class{dataObj}}
#' @keywords classes
#' @examples
#' 
#' showClass("simPopObj")
#' 
#' ## show method: generate an object of class simPop first
#' data(eusilcS)
#' inp <- specifyInput(data=eusilcS, hhid="db030", hhsize="hsize", strata="db040", weight="db090")
#' eusilcP <- simStructure(data=inp, method="direct", basicHHvars=c("age", "rb090"))
#' class(eusilcP)
#' ## shows some basic information:
#' eusilcP
#' 
NULL





#' Population totals Region times Gender for Austria 2006
#' 
#' Population characteristics Region times Gender from Austria.
#' 
#' 
#' @name totalsRG
#' @aliases totalsRG totalsRGtab
#' @docType data
#' @format totalsRG: A data frame with 18 observations on the following 3
#' variables.  \describe{ \item{list("rb090")}{gender; a factor with levels
#' \code{female} \code{male}} \item{list("db040")}{region; a factor with levels
#' \code{Burgenland} \code{Carinthia} \code{Lower Austria,} \code{Salzburg}
#' \code{Styria} \code{Tyrol} \code{Upper Austria} \code{Vienna}
#' \code{Vorarlberg}} \item{list("Freq")}{totals; a numeric vector} }
#' totalsRGtab: a two-dimensional table holding the same information
#' @source StatCube - statistical data base,
#' \url{http://www.statistik.at/web_de/services/datenbank_superstar/}
#' @keywords datasets
#' @examples
#' 
#' data(totalsRG)
#' totalsRG
#' data(totalsRGtab)
#' totalsRGtab
#' 
NULL

setClassUnion('dataframeOrNULL', c('data.frame', 'NULL'))
setClassUnion('numericOrNULL', c('numeric', 'NULL'))
setClassUnion('characterOrNULL', c('character', 'NULL'))
setClassUnion('logicalOrNULL', c('logical', 'NULL'))
setClassUnion('matrixOrNULL', c('matrix', 'NULL'))
setClassUnion('listOrNULL', c('list', 'NULL'))

setClass(
  Class='dataObj',
  representation=representation(
    data='dataframeOrNULL',
    hhid='characterOrNULL',
    pid='characterOrNULL',
    hhsize='characterOrNULL',
    weight='characterOrNULL',
    strata='characterOrNULL',
    ispopulation='logicalOrNULL'
  ),
  prototype=prototype(
    data=NULL,
    hhid=NULL,
    pid=NULL,
    hhsize=NULL,
    weight=NULL,
    strata=NULL,
    ispopulation=NULL
  ),
  validity=function(object) {
    return(TRUE)
  }
)

setClassUnion('dataObjOrNULL', c('dataObj', 'NULL'))

setClass(
  Class='simPopObj',
  representation=representation(
    sample='dataObjOrNULL',
    table='dataframeOrNULL',
    pop='dataObjOrNULL',
    basicHHvars='characterOrNULL'
  ),
  prototype=prototype(
    sample=NULL,
    table=NULL,
    pop=NULL,
    basicHHvars=NULL
  ),
  validity=function(object) {
    return(TRUE)
  }
)


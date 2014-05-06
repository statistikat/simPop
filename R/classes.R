setClassUnion('dataframeOrNULL', c('data.frame', 'NULL'))
setClassUnion('numericOrNULL', c('numeric', 'NULL'))
setClassUnion('characterOrNULL', c('character', 'NULL'))
setClassUnion('logicalOrNULL', c('logical', 'NULL'))
setClassUnion('matrixOrNULL', c('matrix', 'NULL'))
setClassUnion('listOrNULL', c('list', 'NULL'))

setClass(
  Class='popObj',
  representation=representation(
    data='dataframeOrNULL',
    hhid='characterOrNULL',
    hhsize='characterOrNULL',
    pid='characterOrNULL',
    strata='characterOrNULL'
  ),
  prototype=prototype(
    data=NULL,
    hhid=NULL,
    hhsize=NULL,
    pid=NULL,
    strata=NULL
  ),
  validity=function(object) {
    return(TRUE)
  }
)

# extends "popObj" by slot "weight"
setClass(
  Class='sampleObj',
  representation=representation(
    weight='characterOrNULL',
    additional='characterOrNULL'
  ),
  contains="popObj",
  prototype=prototype(
    weight=NULL,
    additional=NULL
  ),
  validity=function(object) {
    return(TRUE)
  }
)

setClassUnion('popObjOrNULL', c('popObj', 'NULL'))
setClassUnion('sampleObjOrNULL', c('sampleObj', 'NULL'))

setClass(
  Class='synthPopObj',
  representation=representation(
    sample='sampleObjOrNULL',
    table='dataframeOrNULL',
    pop='popObjOrNULL'
  ),
  prototype=prototype(
    sample=NULL,
    table=NULL,
    pop=NULL
  ),
  validity=function(object) {
    return(TRUE)
  }
)


# ----------------------------------------
# Authors: Stefan Kraft and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

spMosaic <- function(x, ...) UseMethod("spMosaic")

spMosaic.default <- function(x, weights = NULL, dataS, dataP, ...) {
    if(is.null(x)) tab <- spTable(dataS, dataP, weights=weights)
    else tab <- spTable(dataS, dataP, x, weights)
    spMosaic(tab, ...)
}

## FIXME: does not work yet
#spMosaic.formula <- function(x, weights = NULL, dataS, dataP, ...) {
#    weights <- eval(substitute(weights), dataS, environment(x))  # evaluate
#    tab <- spTable(dataS, dataP, x, weights)
#    spMosaic(tab, ...)
#}

#spMosaic.spTable <- function(x, main=NULL, ...) {
#    # initializations
#    if(is.null(main)) main <- c("Expected", "Realized")
#    # define local version of 'mosaic'
#    localMosaic <- function(x, ..., newpage) {
#        mosaic(x, ..., newpage=FALSE)
#    }
#    # create plot
#    grid.newpage()
#    vp1 <- viewport(layout=grid.layout(nrow=1, ncol=2))  # layout
#    pushViewport(vp1)
#    vp2 <- viewport(layout.pos.row=1, layout.pos.col=1)
#    pushViewport(vp2)
#    localMosaic(x$expected, main=main[1], ...)  # expected (from sample)
#    popViewport()
#    vp2 <- viewport(layout.pos.row=1, layout.pos.col=2)
#    pushViewport(vp2)
#    localMosaic(x$realized, main=main[2], ...)  # realized (population)
#    upViewport()   
#}

spMosaic.spTable <- function(x, ...) {
    # initializations
    tab <- as.table(x)
    dn <- dimnames(tab)
    dn[[length(dn)]] <- c("Sample", "Population")
    names(dn)[length(dn)] <- "Data"
    dimnames(tab) <- dn
    # define local version of 'cotabplot'
    localCotabplot <- function(x, ..., cond, panel) {
        cotabplot(x, cond="Data", ...)
    }
    # produce plot
    localCotabplot(tab, ...)
}

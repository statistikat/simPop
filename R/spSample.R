# ---------------------------------------
# Author: Stefan Kraft
#         Vienna University of Technology
# ---------------------------------------

# 'sample' gives warning for alias sampling due to incompatible results with 
# earlier versions of R

spSample <- function(n, p) {
#    suppressWarnings(sample(length(p), size=n, replace=TRUE, prob=p))
    sample(length(p), size=n, replace=TRUE, prob=p)
}

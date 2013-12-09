# ---------------------------------------
# Author: Stefan Kraft
#         Vienna University of Technology
# ---------------------------------------

spPredict <- function(mod, newdata, intercept = FALSE) {
    #beta <- mod$coef
    beta <- coef(mod)
    if(intercept) X <- cbind(1 , newdata) else X <- newdata 
    as.numeric(X %*% beta)
}

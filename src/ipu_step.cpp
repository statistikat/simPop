#include <Rcpp.h>
using namespace Rcpp;

//' Perform one step of iterative proportional updating
//' 
//' C++ routines to invoke a single iteration of the Iterative proportional updating (IPU) scheme. Targets and classes 
//' are assumed to be one dimensional in the `ipu_step` functions. `combine_factors` aggregates several vectors of 
//' type factor into a single one to allow multidimensional ipu-steps. See examples.
//' 
//' `ipu_step` returns the adjusted weights. `ipu_step_ref` does the same, but updates `w` by reference rather than 
//' returning. `ipu_step_f` returns a multiplicator: adjusted weights divided by unadjusted weights. `combine_factors` is
//' designed to make `ipu_step` work with contingency tables produced by [xtabs].
//' 
//' @md
//' @name ipu_step
//' @param w An numeric vector of weights. All entries should be positive.
//' @param classes A factor variable. Must have the same length as `w`.
//' @param targets key figure to target with the ipu scheme. A numeric verctor of the same length as `levels(classes)`. 
//' This can also be a `table` produced by `xtabs`. See examples.
//' @examples
//' 
//' ############# one-dimensional ipu ##############
//' 
//' ## create random data
//' nobs <- 10
//' classLabels <- letters[1:3]
//' dat = data.frame(
//'   weight = exp(rnorm(nobs)),
//'   household = factor(sample(classLabels, nobs, replace = TRUE))
//' )
//' dat
//' 
//' ## create targets (same lenght as classLabels!)
//' targets <- 3:5
//' 
//' ## calculate weights
//' new_weight <- ipu_step(dat$weight, dat$household, targets)
//' cbind(dat, new_weight)
//' 
//' ## check solution
//' xtabs(new_weight ~ dat$household)
//' 
//' ## calculate weights "by reference"
//' ipu_step_ref(dat$weight, dat$household, targets)
//' dat
//' 
//' ############# multidimensional ipu ##############
//' 
//' ## load data
//' factors <- c("time", "sex", "smoker", "day")
//' data(tips, package = "reshape2")
//' tips <- tips[factors]
//' 
//' ## combine factors
//' cf <- combine_factors(tips, names(tips))
//' cbind(tips, cf)[sample(nrow(tips), 10),]
//' 
//' ## adjust weights
//' con <- xtabs(~., tips)
//' weight <- rnorm(nrow(tips)) + 5
//' adjusted_weight <- ipu_step(weight, cf, con)
//' 
//' ## check outputs
//' con2 <- xtabs(adjusted_weight ~ ., data = tips)
//' sum((con - con2)^2)
//' 
//' @rdname ipu_step
//' @export
// [[Rcpp::export]]
void ipu_step_ref(NumericVector w, IntegerVector classes, NumericVector targets) {
  CharacterVector levels(classes.attr("levels"));
  int nclasses = levels.size();
  if(targets.length() != nclasses)
    stop("number of levels does not match the length of targets");
  NumericVector targets2(targets.length());
  for(int i = 0; i < w.size(); i++){
    int cl = classes[i] - 1;
    targets2[cl] += w[i];
  }
  for(int i = 0; i < w.size(); i++){
    int cl = classes[i] - 1;
    w[i] *= targets[cl];
    w[i] /= targets2[cl];
  }
}

//' @rdname ipu_step
//' @export
// [[Rcpp::export]]
NumericVector ipu_step(NumericVector w, IntegerVector classes, NumericVector targets){
  NumericVector w_copy( Rcpp::clone( w ) );
  ipu_step_ref(w_copy, classes, targets);
  return w_copy;
}

//' @rdname ipu_step
//' @export
// [[Rcpp::export]]
NumericVector ipu_step_f(NumericVector w, IntegerVector classes, NumericVector targets){
  CharacterVector levels(classes.attr("levels"));
  int nclasses = levels.size();
  if(targets.length() != nclasses)
    stop("number of levels does not match the length of targets");
  NumericVector targets2(targets.length());
  for(int i = 0; i < w.size(); i++){
    int cl = classes[i] - 1;
    targets2[cl] += w[i];
  }
  NumericVector adjustments(classes.size());
  for(int i = 0; i < classes.size(); i++){
    int cl = classes[i] - 1;
    adjustments[i] = targets[cl]/targets2[cl];
  }
  return adjustments;
}
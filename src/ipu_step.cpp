#include <Rcpp.h>
using namespace Rcpp;

//' Perform one step of iterative proportional updating
//' 
//' C++ routine to invoke a single iteration of the IPU scheme. Targets and classes are assumed to be one dimensional. 
//' `ipu_step` returns the adjusted weigth verctor while `ipu_step_ref` updates the input `w` by reference.
//' 
//' @md
//' @name ipu_step
//' @param w An numeric vector of weights. All entries should be positive.
//' @param classes A factor variable. Must have the same length as `w`.
//' @param targets key figure to target with the ipu scheme. A numeric verctor of the same length as `levels(classes)`.
//' @examples
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
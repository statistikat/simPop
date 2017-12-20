#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' Calculate mean by factors
//' 
//' These functions calculate the arithmetic and geometric mean of the weight for each class. `geometric_mean` and
//' `arithmetic_mean` return a `numeric` vector of the same length as `w` which stores the averaged weight for each 
//' observation. `geometric_mean_reference` returns the same value by reference, i.e. the input value `w` gets 
//' overwritten by the updated weights. See examples.
//' 
//' @md
//' @name cpp_mean
//' @param w An numeric vector. All entries should be positive.
//' @param classes A factor variable. Must have the same length as `w`.
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
//' ## calculate weights with geometric_mean
//' geom_weight <- geometric_mean(dat$weight, dat$household)
//' cbind(dat, geom_weight)
//' 
//' ## calculate weights with arithmetic_mean
//' arith_weight <- arithmetic_mean(dat$weight, dat$household)
//' cbind(dat, arith_weight)
//' 
//' ## calculate weights "by reference"
//' geometric_mean_reference(dat$weight, dat$household)
//' dat
//' 
//' @rdname cpp_mean
//' @export
// [[Rcpp::export]]
void geometric_mean_reference(NumericVector& w, const IntegerVector& classes) {
  CharacterVector levels(classes.attr("levels"));
  int nclasses = levels.size();
  
  NumericVector sums( nclasses );
  NumericVector sizes( nclasses );
  
  for (int i = 0; i < w.size(); i++){
    int cl = classes[i] - 1;
    sums[cl] += log(w[i]);
    sizes[cl] += 1;
  }
  
  NumericVector means(nclasses);
  for(int cl = 0; cl < means.size(); cl++){
    means[cl] = exp(sums[cl]/sizes[cl]);
  }
  
  for(int i = 0; i < w.size(); i++){
    int cl = classes[i] - 1;
    w[i] = means[cl];
  }
}
//' @rdname cpp_mean
//' @export
// [[Rcpp::export]]
NumericVector geometric_mean( const NumericVector& w, const IntegerVector& classes ){
  NumericVector w_copy( Rcpp::clone( w ) );
  geometric_mean_reference( w_copy, classes );
  return w_copy;
}

//' @rdname cpp_mean
//' @export
// [[Rcpp::export]]
NumericVector arithmetic_mean(const NumericVector& w, const IntegerVector& classes) {
  CharacterVector levels(classes.attr("levels"));
  int nclasses = levels.size();
  
  NumericVector sums( nclasses );
  NumericVector sizes( nclasses );
  NumericVector w_copy(Rcpp::clone(w));
  
  for (int i = 0; i < w.size(); i++){
    int cl = classes[i] - 1;
    sums[cl] += w[i];
    sizes[cl] += 1;
  }
  
  NumericVector means(nclasses);
  for(int cl = 0; cl < means.size(); cl++){
    means[cl] = sums[cl]/sizes[cl];
  }
  
  for(int i = 0; i < w.size(); i++){
    int cl = classes[i] - 1;
    w_copy[i] = means[cl];
  }
  return w_copy;
}


#include <Rcpp.h>
using namespace Rcpp;

//' @rdname computeFrac
//' @export
// [[Rcpp::export]]
NumericVector computeLinear(double curValue, 
                            double target, 
                            const NumericVector& x,
                            const NumericVector& w, 
                            double boundLinear = 10) {
  double h = 0.0;
  double j = 0.0;
  double N = 0.0;
  
  for(int i = 0; i < x.size(); i++){
    h += w[i]*x[i];
    j += w[i]*x[i]*x[i];
    N += w[i];
  }
  
  double b = (target-N*j/h)/(h-N*j/h);
  double a = (N-b*N)/h;

  NumericVector f(x.size());
  for(int i = 0; i < x.size(); i++)
    f[i] = a*x[i] + b;
  
  //apply bounds
  for(int i = 0; i < x.size(); i++){
    if (f[i] < 1.0/boundLinear)
      f[i] = 1.0/boundLinear;
    if (f[i] > boundLinear)
      f[i] = boundLinear;
  }
  
  return f;
}

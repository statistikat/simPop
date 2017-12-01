#include <Rcpp.h>
using namespace Rcpp;

//' @rdname ipu2
//' @export
// [[Rcpp::export]]
NumericVector computeLinear(double curValue, 
                            double Value, 
                            const NumericVector& numericVar,
                            const NumericVector& weightVec, 
                            double boundLinear = 10) {
  double h = 0.0;
  double j = 0.0;
  double N = 0.0;
  
  for(int i = 0; i < numericVar.size(); i++){
    h += weightVec[i]*numericVar[i];
    j += weightVec[i]*numericVar[i]*numericVar[i];
    N += weightVec[i];
  }
  
  double b = (Value-N*j/h)/(h-N*j/h);
  double a = (N-b*N)/h;

  NumericVector f(numericVar.size());
  for(int i = 0; i < numericVar.size(); i++)
    f[i] = a*numericVar[i] + b;
  
  //apply bounds
  for(int i = 0; i < numericVar.size(); i++){
    if (f[i] < 1.0/boundLinear)
      f[i] = 1.0/boundLinear;
    if (f[i] > boundLinear)
      f[i] = boundLinear;
  }
  
  return f;
}

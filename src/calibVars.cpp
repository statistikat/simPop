#include <R.h>
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerMatrix binary_representation(IntegerVector levels, IntegerVector values) {
  IntegerMatrix out(values.size(), levels.size());
  for ( int j=0; j<levels.size(); ++j ) {
    for ( int i=0; i<values.size(); ++i ) {
      if ( levels[j] == values[i] ) {
        out(i,j) = 1;
      }
    }
  }
  return(out);
}

#include <R.h>
#include <Rcpp.h>
using namespace Rcpp;

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

IntegerMatrix binary_representation(IntegerVector levels, IntegerVector values);
RcppExport SEXP synthPop_binary_representation(SEXP levelsSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type levels(levelsSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type values(valuesSEXP );
        IntegerMatrix __result = binary_representation(levels, values);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
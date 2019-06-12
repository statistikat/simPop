#include <R.h>
#include <cmath>
#include <Rcpp.h>

void R_init_simPop(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector ipu_work(NumericMatrix inp, NumericVector con, NumericVector w, double eps, IntegerVector verbose) {
  int nr_con = con.size();
  int nr_rows = inp.nrow();
  NumericVector gamma_vals(nr_con);
  NumericVector gamma_vals_new(nr_con);
  double gamma, gamma_new, delta, v;
  bool verb = false;
  if ( verbose[0] == 1 ) {
    verb = true;
  }

  for ( int i=0; i < nr_con; ++i ) {
    gamma_vals[i] = (fabs(sum(inp(_,i)*w)-con[i])) / con[i];
  }
  gamma = mean(gamma_vals);

  bool run_ind = true;
  int counter = 0;
  while ( run_ind ) {
    counter = counter + 1;
    for ( int j=0; j < nr_con; ++j ) {
      v = con[j] / sum(inp(_,j)*w);
      for ( int k=0; k<nr_rows; ++k ) {
        if ( inp(k,j) != 0 ) {
          w[k] = v*w[k];
        }
      }
    }
    //double meanweight;
    //meanweight=mean(w);
    // recalculate gamma_vals
    for ( int i=0; i < nr_con; ++i ) {
      gamma_vals_new[i] = (fabs(sum(inp(_,i)*w)-con[i])) / con[i];
    }
    gamma_new = mean(gamma_vals_new);

    delta = fabs(gamma_new - gamma);
    if ( verb ) {
      Rprintf("improvement in run %d: %g | gamma_new=%g | gamma=%g \n", counter, delta, gamma_new,gamma);
    }

    if ( gamma_new < eps ) {
      if ( verb ) {
        Rprintf("ipu finished after %d interations!\n", counter);
      }
      run_ind = false;
    } else if ( delta < eps/10 ) {
      if ( verb ) {
        Rprintf("WARNING: not converted \n");
      }
        run_ind = false;
    }else {
      for ( int k=0; k<nr_con; ++k ) {
        gamma_vals[k] = gamma_vals_new[k];
      }
      gamma = gamma_new;
    }
  }
  return(w);
}

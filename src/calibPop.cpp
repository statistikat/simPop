#include <R.h>
#include <cmath>
#include <cstdlib>
#include <Rcpp.h>
using namespace Rcpp;

// objective function
double calc_obj(IntegerMatrix inp, IntegerVector weights, NumericVector totals) {
  // in jeder zeile von inp stehen die indices eines constraints
  // der dazugehoerige totalwert steht in totals[i]
  // weights == 1 -> person ausgewaehlt  
  double obj=0.0;
  for ( int i=0; i < totals.size(); ++i ) {
    obj = obj + fabs(sum(inp(i,_)*weights) - totals[i]);
  }
  return(obj);  
}

// generate a new random solution by setting some households from 0->1 and some from1->0
IntegerVector generate_new_solution(IntegerVector weights, IntegerVector hh_ids, 
  IntegerVector hh_head, double factor) {
  std::vector<int> hh_ids_unique;
  std::vector<int> hh_status_unique;
  for ( int i=0; i<hh_ids.size(); ++i ) {
    if ( hh_head[i] == 1 ) {
      hh_ids_unique.push_back(hh_ids[i]);
      hh_status_unique.push_back(weights[i]);
    }
  }
  int nr_draws = floor(factor);
  
  int nr_hh = hh_status_unique.size();
  int nr_active = std::accumulate(hh_status_unique.begin(), hh_status_unique.end(), 0);
  int nr_inactive = nr_hh - nr_active;
  
  IntegerVector indices_active(nr_active);
  IntegerVector indices_inactive(nr_active);
  
  int counter_act = 0;
  int counter_inact = 0;
  for ( int i=0; i<nr_hh; ++i ) {
    // active
    if ( hh_status_unique[i] == 1 ) {      
      indices_active[counter_act] = hh_ids_unique[i];
      counter_act = counter_act + 1;
    }
    // inactive
    if ( hh_status_unique[i] == 0 ) {
      indices_inactive[counter_inact] = hh_ids_unique[i];
      counter_inact = counter_inact + 1;
    }
  }
  
  // sampling
  IntegerVector draws1(nr_draws);
  IntegerVector draws2(nr_draws);
  
  // fixme: sample without replacement
  int s = 0;
  for ( int i=0; i < nr_draws; ++i ) {
    // select an active household and set it to inactive
    s = rand() % nr_active;     // v2 in the range 0 to nr_active
    draws1[i] = indices_active[s];
        
    // select an inactive household and set it to active
    s = rand() % nr_inactive;   // v2 in the range 0 to nr_inactive
    draws2[i] = indices_inactive[s];    
  }
  
  // fixme: std::find? | rcpp match?
  IntegerVector w_neu(weights.size());
  for ( int i=0; i<nr_draws; ++i ) {
    for ( int k=0; k<weights.size(); ++k ) {
      w_neu[k] = weights[k];
      if ( hh_ids[k] == draws1[i] ) {
        w_neu[k] = 0;
      } else if ( hh_ids[k] == draws2[i] ) {
        w_neu[k] = 1;
      }      
    }
  }  
  return(w_neu);
}

double median_rcpp(IntegerVector x) {
   IntegerVector y = clone(x);
   int n, half;
   double y1, y2;
   n = y.size();
   half = n / 2;
   if(n % 2 == 1) {
      // median for odd length vector
      std::nth_element(y.begin(), y.begin()+half, y.end());
      return y[half] / 1.0;
   } else {
      // median for even length vector
      std::nth_element(y.begin(), y.begin()+half, y.end());
      y1 = y[half];
      std::nth_element(y.begin(), y.begin()+half-1, y.begin()+half);
      y2 = y[half-1];
      return (y1 + y2) / 2.0;
   }
}

double calc_factor(double obj, IntegerVector hh_head, IntegerVector hh_size, IntegerVector weights) {
  std::vector<int> hh_size_unique; // unique and active
  for ( int i=0; i<weights.size(); ++i ) {
    if ( hh_head[i] == 1 && weights[i] == 1 ) {
      hh_size_unique.push_back(hh_size[i]);
    }
  }  
  
  int nr_hh = hh_size_unique.size();
  IntegerVector sizes(nr_hh);
  for ( int i=0; i<nr_hh; ++i ) {
    sizes[i] = hh_size_unique[i];
  }
  double factor = obj / median_rcpp(sizes)/10;
  return(factor);
}  

// [[Rcpp::export]]
IntegerVector calibPop_work(IntegerMatrix inp, NumericVector totals, IntegerVector weights, 
  List hh_info, List params) {
    
  IntegerVector new_weights(weights.size());
  
  int nr_con = totals.size();  
  bool verb = false;
  if ( params["verbose"] == 1 ) {
    verb = true;
  }

  IntegerVector hh_ids = hh_info["hh_ids"];
  IntegerVector hh_head = hh_info["hh_head"];
  IntegerVector hh_size = hh_info["hh_size"];

  double temp = params["temp"];
  double eps_factor = params["eps_factor"];
  double maxiter = params["maxiter"];
  double temp_cooldown = params["temp_cooldown"];
  double factor_cooldown = params["factor_cooldown"];
  double min_temp = params["min_temp"];
  int cooldown = 0; 
  double obj, obj_new = 0.0;
  
  //acceptable error
  double eps = eps_factor*sum(weights);
    
  // calculate objective value based on initial solution
  obj = calc_obj(inp, weights, totals);
  
  double factor = calc_factor(obj, hh_head, hh_size, weights);
  if ( obj <= eps ) {
    Rprintf("nothing to do, already finished! obj=%g | eps=%g\n", obj, eps);
    return(weights);
  }
  
  int counter = 0;
  int change = 0;
  double prob = 0.0;
  NumericVector samp_result(1);
  while( temp > min_temp ) {  
    //Rprintf("temp=%f | min_temp=%f\n", temp, min_temp);
    counter = 1;
    while( counter < maxiter ) {
      // swap zeros to ones and ones to zeros
      new_weights = generate_new_solution(weights, hh_ids, hh_head, factor);
      
      // calculate new objective value based on new solution
      obj_new = calc_obj(inp, new_weights, totals);
      if ( verb ) {
        Rprintf("obj_new=%g | obj_old=%g\n", obj_new, obj);
      }
      
      if ( obj_new <= eps ) {
        obj = obj_new;
        for ( int z=0; z<weights.size(); ++z ) {
          weights[z] = new_weights[z];
        }
        change = change+1;
        Rprintf("breaking: obj_new: %g | eps=%g\n", obj_new, eps);
        break;
      }      
      // if new solution is better than old one, we accept
      if ( obj_new <= obj ) { 
        for ( int z=0; z<weights.size(); ++z ) {
          weights[z] = new_weights[z];
        }
        obj = obj_new;
        change = change + 1;
      } 
      // if new solution is worse than current, we accept based on temp and value
      if ( obj_new > obj ) {      
        prob = exp(-(obj_new-obj)/temp);
        samp_result = runif(1);
        if ( samp_result[0] <= prob ) { 
          for ( int z=0; z<weights.size(); ++z ) {
            weights[z] = new_weights[z];
          }
          obj = obj_new;
          change = change + 1;
        }
      }
      counter = counter + 1; 
    }
    // decrease temp and decrease factor accordingly
    // decrease temp by a const fraction (simple method used for testing only)
    temp = temp_cooldown*temp;
    factor = floor(factor_cooldown*factor);
    if ( factor == 0 ) {
      factor = 1;
    }
    cooldown = cooldown + 1;
    if ( obj_new <= eps | cooldown == 500 ) {
      Rprintf("breaking: obj_new: %g | eps=%g\n", obj_new, eps);
      break;
    }    
  }
  return(weights);
}

// exporting for package (compileAttributes())
IntegerVector calibPop_work(IntegerMatrix inp, NumericVector totals, IntegerVector weights, List hh_info, List params);
RcppExport SEXP synthPop_calibPop_work(SEXP inpSEXP, SEXP totalsSEXP, SEXP weightsSEXP, SEXP hh_infoSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type inp(inpSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type totals(totalsSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type weights(weightsSEXP );
        Rcpp::traits::input_parameter< List >::type hh_info(hh_infoSEXP );
        Rcpp::traits::input_parameter< List >::type params(paramsSEXP );
        IntegerVector __result = calibPop_work(inp, totals, weights, hh_info, params);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

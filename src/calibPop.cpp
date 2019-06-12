#include <R.h>
#include <cmath>
#include <cstdlib>
#include <Rcpp.h>
using namespace Rcpp;

// objective function
NumericVector calc_obj(IntegerMatrix inp, IntegerVector weights, NumericVector totals) {
  // in each row of inp we have the indices of a single constraint
  // the corresponding total is available from totals[i]
  // weights == 1 -> person is in synthetic population
  NumericVector obj(1);
  for ( int i=0; i < totals.size(); ++i ) {
    obj[0] = obj[0] + fabs(sum(inp(i,_)*weights) - totals[i]);
  }
  return(obj);
}

// generate a new random solution by setting some households from 0->1 and some from 1->0
IntegerVector generate_new_solution(IntegerVector weights, IntegerVector hh_ids, IntegerVector hh_head, double factor) {
  std::vector<int> hh_ids_unique;
  std::vector<int> hh_status_unique;
  for ( int i=0; i<hh_ids.size(); ++i ) {
    if ( hh_head[i] == 1 ) {
      hh_ids_unique.push_back(hh_ids[i]);
      hh_status_unique.push_back(weights[i]);
    }
  }

  int nr_draws = (int)floor(factor);
  int nr_hh = hh_status_unique.size();
  int nr_active = std::accumulate(hh_status_unique.begin(), hh_status_unique.end(), 0);
  int nr_inactive = nr_hh - nr_active;

  std::vector<int> indices_active, indices_inactive;
  for ( int i=0; i<nr_hh; ++i ) {
    if ( hh_status_unique[i] == 1 ) {
      // active
      indices_active.push_back(hh_ids_unique[i]);
    } else {
      // inactive
      indices_inactive.push_back(hh_ids_unique[i]);
    }
  }

  // sampling
  std::vector<int> draws1, draws2;

  // fixme: sample without replacement
  int s = 0;
  NumericVector ss(1);
  for ( int i=0; i < nr_draws; ++i ) {
    // select an active household and set it to inactive
    ss = floor(runif(1,0,nr_active));     // range: 0 to nr_active
    s = (int)ss[0];
    if ( s > nr_active ) {
      s = nr_active;
    }
    draws1.push_back(indices_active[s]);

    // select an inactive household and set it to active
    ss = floor(runif(1,0,nr_inactive));
    s = (int)ss[0];
    if ( s > nr_inactive ) {
      s = nr_inactive;
    }
    draws2.push_back(indices_inactive[s]);
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

//[[Rcpp::export]]
IntegerVector calibPop_work(IntegerMatrix inp, NumericVector totals, IntegerVector weights,
  List hh_info, List params) {

  IntegerVector new_weights(weights.size());

  bool verb = false;
  int verbTF = params["verbose"];
  if ( verbTF == 1 ) {
    verb = true;
  }

  IntegerVector hh_ids = hh_info["hh_ids"];
  IntegerVector hh_head = hh_info["hh_head"];
  IntegerVector hh_size = hh_info["hh_size"];
  NumericVector median_hhsize_x = hh_info["median_hhsize"];
  double median_hhsize = median_hhsize_x[0];
  NumericVector temp_x = params["temp"];
  double temp = temp_x[0];
  NumericVector eps_factor_x = params["eps_factor"];
  double eps_factor = eps_factor_x[0];
  NumericVector maxiter_x = params["maxiter"];
  double maxiter = maxiter_x[0];
  NumericVector temp_cooldown_x = params["temp_cooldown"];
  double temp_cooldown = temp_cooldown_x[0];
  NumericVector factor_cooldown_x = params["factor_cooldown"];
  double factor_cooldown = factor_cooldown_x[0];
  NumericVector min_temp_x = params["min_temp"];
  double min_temp = min_temp_x[0];
  int cooldown = 0;
  double obj,obj_new = 0.0;

  //acceptable error
  double eps = eps_factor*sum(weights);

  // calculate objective value based on initial solution
  obj = calc_obj(inp, weights, totals)[0];

  double factor = (median_hhsize / 5) * obj;
  if ( obj <= eps ) {
    if ( verb ) {
      Rprintf("We have nothing to do and are already finished!\nValue of objective function: %g | (required precision=%g)\n", obj, eps);
    }
    return(weights);
  }

  int counter = 0;
  int change = 0;
  double prob = 0.0;
  NumericVector samp_result(1);
  while ( temp > min_temp ) {
    if ( verb ) {
      Rprintf("current temperature: %f (minimal temp=%f)\n", temp, min_temp);
    }
    counter = 1;
    while( counter < maxiter ) {
      // swap zeros to ones and ones to zeros
      new_weights = generate_new_solution(weights, hh_ids, hh_head, factor);

      // calculate new objective value based on new solution
      obj_new = calc_obj(inp, new_weights, totals)[0];

      if ( verb ) {
        Rprintf("--> value of objective: %g (previous=%g)\n", obj_new, obj);
      }
      if ( obj_new <= eps ) {
        obj = obj_new;
        for ( int z=0; z<weights.size(); ++z ) {
          weights[z] = new_weights[z];
        }
        change = change+1;
        if ( verb ) {
          Rprintf("Required precision reached!\nValue of objective function: %g (required precision=%g)\n", obj_new, eps);
        }
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
          if ( verb ) {
            Rprintf("----> worse solution accepted! (samp=%g | prob=%g)]\n", samp_result[0], prob);
          }
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
    if ( (obj_new <= eps) | (cooldown == 500) ) {
      if ( verb ) {
        Rprintf("Required precision reached!\nValue of objective function: %g (required precision=%g)\n", obj_new, eps);
      }
      break;
    }
  }
  return(weights);
}

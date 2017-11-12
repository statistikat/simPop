#include <Rcpp.h>
using namespace Rcpp;
// Author: Kirill Müller https://github.com/krlmlr/wrswoR
#include <queue>

void check_args(int n, int size, const NumericVector& prob) {
  if (n < size) {
    Rcpp::stop("cannot take a sample larger than the population");
  }
  
  if (prob.length() != n) {
    Rcpp::stop("incorrect number of probabilities");
  }
}


template <class T>
T _rexp_divide_by(T t) { return Rf_rexp(1.0) / t; }

struct Indirection{
  Indirection(const Rcpp::IntegerVector& v ) : _v(v) {}
  // Inverted comparison!
  int operator ()(int a) { return _v[a]; }
  const Rcpp::IntegerVector& _v;
};

// [[Rcpp::export(sample_int_expj)]]
IntegerVector sample_int_expj(int n, int size, NumericVector prob) {
  check_args(n, size, prob);
  
  // Corner case
  if (size == 0)
    return IntegerVector();
  
  // Step 1: The first m items of V are inserted into R
  // Step 2: For each item v_i ∈ R: Calculate a key k_i = u_i^(1/w),
  // where u_i = random(0, 1)
  // (Modification: Calculate and store -log k_i = e_i / w where e_i = exp(1),
  //  reservoir is a priority queue that pops the *maximum* elements)
  std::priority_queue<std::pair<double, int> > R;
  
  for (NumericVector::iterator iprob = prob.begin();
       iprob != prob.begin() + size; ++iprob) {
    double k_i = _rexp_divide_by<double>(*iprob);
    R.push(std::make_pair(k_i, iprob - prob.begin() + 1));
  }
  
  // Step 4: Repeat Steps 5–10 until the population is exhausted
  {
    // Step 3: The threshold T_w is the minimum key of R
    // (Modification: This is now the logarithm)
    // Step 10: The new threshold T w is the new minimum key of R
    const std::pair<double, int>& T_w = R.top();
    
    // Incrementing iprob is part of Step 7
    for (NumericVector::iterator iprob = prob.begin() + size; iprob != prob.end(); ++iprob) {
      
      // Step 5: Let r = random(0, 1) and X_w = log(r) / log(T_w)
      // (Modification: Use e = -exp(1) instead of log(r))
      double X_w = Rf_rexp(1.0) / T_w.first;
      
      // Step 6: From the current item v_c skip items until item v_i, such that:
      double w = 0.0;
      
      // Step 7: w_c + w_{c+1} + ··· + w_{i−1} < X_w <= w_c + w_{c+1} + ··· + w_{i−1} + w_i
      for (; iprob != prob.end(); ++iprob) {
        w += *iprob;
        if (X_w <= w)
          break;
      }
      
      // Step 7: No such item, terminate
      if (iprob == prob.end())
        break;
      
      // Step 9: Let t_w = T_w^{w_i}, r_2 = random(t_w, 1) and v_i’s key: k_i = (r_2)^{1/w_i}
      // (Mod: Let t_w = log(T_w) * {w_i}, e_2 = log(random(e^{t_w}, 1)) and v_i’s key: k_i = -e_2 / w_i)
      double t_w = -T_w.first * *iprob;
      double e_2 = std::log(Rf_runif(std::exp(t_w), 1.0));
      double k_i = -e_2 / *iprob;
      
      // Step 8: The item in R with the minimum key is replaced by item v_i
      R.pop();
      R.push(std::make_pair(k_i, iprob - prob.begin() + 1));
    }
  }
  
  IntegerVector ret(size);
  
  for (IntegerVector::iterator iret = ret.end(); iret != ret.begin(); ) {
    --iret;
    
    if (R.empty()) {
      stop("Reservoir empty before all elements have been filled");
    }
    
    *iret = R.top().second;
    R.pop();
  }
  
  if (!R.empty()) {
    stop("Reservoir not empty after all elements have been filled");
  }
  
  return ret;
}
#include <Rcpp.h>
using namespace Rcpp;
//' @rdname ipu_step
//' @name ipu_step
//' @md
//' 
//' @param dat a data.frame containing the factor variables to be combined.
//' @param factor_columns a `character` vector containing the column names in `dat` to be combined. All columns must
//' be factors.
//' @param as_factor Wheter to return an integer or a factor variable.
//' 
//' @export
// [[Rcpp::export]]
IntegerVector combine_factors(DataFrame& dat, CharacterVector factor_columns, bool asfactor = true) {
  IntegerVector new_column(dat.rows());
  int nlevelsout = 1;
  
  // calculate indices and number of possible values
  for(int j = 0; j < factor_columns.size(); j++){
    String column_name = factor_columns[factor_columns.size() - j - 1];
    IntegerVector column = dat[column_name];
    CharacterVector levels = column.attr("levels");
    int nlevels = levels.size();
    
    for(int i = 0; i < dat.rows(); i++){
      new_column[i] *= nlevels;
      new_column[i] += column[i] - 1;
    }
    
    nlevelsout *= nlevels;
  }
  
  // correct for R indexing (necessary for factor variables)
  for(int i = 0; i < dat.rows(); i++)
    new_column[i] += 1;
  
  if(asfactor){
    IntegerVector levs(nlevelsout);
    for(int i = 0; i < nlevelsout; i++)
      levs[i] = i + 1;
    new_column.attr("levels") = as<CharacterVector>(levs);
    new_column.attr("class") = "factor";
  }
  
  // return
  return new_column;
}

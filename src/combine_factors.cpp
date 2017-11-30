#include <Rcpp.h>
using namespace Rcpp;

//' Combine multiple factorvariables to a single one
//' 
//' In order to apply [ipu_step] to multidimensional contigency-tables, this function is used to reduce
//' several factor variables to a single factor variable. The number of levels for the result is the 
//' product of the number of levels of the input factors.
//' 
//' @md
//' @param dat A `data.frame` object that holds several factor variables.
//' @param factor_columns A `character` naming all factors to be combined
//' @param asfactor A `logical` of length one, specifying whether the return value should be a `factor` (default)
//' or an `intger`.
//' @examples
//' 
//' ################## create data ##################
//' factors <- c("time", "sex", "smoker", "day")
//' tips <- reshape2::tips[factors]
//' 
//' ################ combine factors ################
//' cf <- combine_factors(tips, names(tips))
//' cbind(tips, cf)[sample(nrow(tips), 10),]
//' 
//' ############# high dimensional ipu ##############
//' con <- xtabs(~., tips)
//' weight <- rnorm(nrow(tips)) + 5
//' adjusted_weight <- ipu_step(weight, cf, con)
//' 
//' #################### check ######################
//' con2 <- xtabs(adjusted_weight ~ ., data = tips)
//' sum((con - con2)^2)
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

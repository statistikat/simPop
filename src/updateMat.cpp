#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix updateMatC(NumericMatrix M,NumericVector add_row, NumericVector add_col, NumericVector remove_row, NumericVector remove_col, NumericVector hhsize, NumericVector hhid) {
  
  // define Variables
  int N = M.nrow();
  int n = add_col.size();
  int size=0;
  int id=0;
  int j_low=0; //for internal loop
  int j_high=0; // for internal loop
  int count_size=0; //for internal loop
  int col_i=0; // for internal loop
  
  /* 
  not working yet(!!!)
  NumericVector remove_row(n);
  NumericVector add_row(n);
  
  // transform index_add values and index_remove values
  // define row for removal and adding
  for(int i=0;i++;i<n){
    remove_row[i] = index_remove[i]%ncol -1;
    add_row[i] = index_add[i]%ncol-1;
  }
  remove_row[remove_row==-1] = nrow;
  add_row[add_row==-1] = nrow;
  
  // define cols for removal and adding
  NumericVector remove_col = floor(index_remove/(nrow+1));
  NumericVector add_col = floor(index_add/(nrow+1));
  
  size = hhsize[remove_row[0]];
  id = hhid[remove_row[0]];
  
  j_low = remove_row[0]-(size-1);
  j_low = std::max(0,j_low);
  j_high = remove_row[0]+(size-1);
  j_high = std::min(j_high,N-1);
  
  count_size =0; //set counter to zero
  NumericVector test(N);
  col_i = remove_col[0];
  for(int j=j_low;j<(j_high+1);j++){
    if(hhid[j]==id){
      M(j,col_i) = 0;
      count_size = count_size + 1;
    }
    test[j] = hhid[j];
    if(count_size==size){
      break;
    }
  }
  //return test;
  
  return M;
  */
  
  // go through Matrix to add and remove households
  for(int i=0;i<n;i++){
    // remove household
    size = hhsize[remove_row[i]];
    id = hhid[remove_row[i]];
    
    j_low = remove_row[i]-(size-1);
    j_low = std::max(0,j_low);
    j_high = remove_row[i]+(size-1);
    j_high = std::min(j_high,N-1);
    
    count_size =0; //set counter to zero
    col_i = remove_col[i];
    for(int j=j_low;j<(j_high+1);j++){
      if(hhid[j]==id){
        M(j,col_i) = 0;
        count_size = count_size + 1;
      }
      if(count_size==size){
        break;
      }
    }
    
    // add household
    size = hhsize[add_row[i]];
    id = hhid[add_row[i]];
    
    j_low = add_row[i]-(size-1);
    j_low = std::max(0,j_low);
    j_high = add_row[i]+(size-1);
    j_high = std::min(j_high,N-1);
    
    count_size =0; //set counter to zero
    col_i = add_col[i];
    for(int j=j_low;j<(j_high+1);j++){
      if(hhid[j]==id){
        M(j,col_i) = 1;
        count_size = count_size + 1;
      }
      if(count_size==size){
        break;
      }
    }
  }
  
  // return updated matrix
  return M;
}

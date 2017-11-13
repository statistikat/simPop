#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export(updateVecC)]]
IntegerVector updateVecC(IntegerVector init_weight,IntegerVector add_index, IntegerVector remove_index, IntegerVector hhsize, IntegerVector hhid,int sizefactor) {
  
  // define Variables
  int n = hhsize.size();
  int n_remove = remove_index.size();
  int n_add = add_index.size();
  int n_ar = std::max(n_remove,n_add);
  int size=0;
  int id=0;
  int j_low=0; //for internal loop
  int j_high=0; // for internal loop
  int count_size=0; //for internal loop
  int id_index=0; //for internal loop

  // transform
  IntegerVector add_id(n_add);
  IntegerVector remove_id(n_remove);
  
  // transform index_add values and index_remove values
  // define row for removal and adding
  for(int i=0;i<n_ar;i++){
    if(i<n_remove){
      remove_id[i] = remove_index[i] % n;
    }
    if(i<n_add){
      add_id[i] = add_index[i] % n;
    }
  }
  //return add_id;
  // got through add_index and remove_index
  // and select update init_weight using hhsize and hhid
  for(int i=0;i<n_ar;i++){
    
    if(i<n_remove){
      //remove household 
      // select size and id using modulus operator
      size = hhsize[remove_id[i]];
      id = hhid[remove_id[i]];
      
      j_low = remove_index[i]-(size-1);
      j_low = std::max(remove_index[i]-remove_index[i]%n,j_low);
      j_high = remove_index[i]+(size-1);
      j_high = std::min(j_high,remove_index[i]+n-remove_index[i]%n);
      
      count_size =0; //set counter to zero
      
      for(int j=j_low;j<(j_high+1);j++){
        
        id_index = j % n;
        
        // replace init_weight only if hhid at id_index has same household id
        if(hhid[id_index]==id){
          init_weight[j] = 0;
          count_size = count_size + 1;
        }
        if(count_size==size){
          break;
        }
      }
    }
    
    if(i<n_add){
      // add household
      // select size and id using modulus operator
      size = hhsize[add_id[i]];
      id = hhid[add_id[i]];
      
      j_low = add_index[i]-(size-1);
      j_low = std::max(add_index[i]-add_index[i]%n,j_low);
      j_high = add_index[i]+(size-1);
      j_high = std::min(j_high,add_index[i]+n-add_index[i]%n);
      
      count_size =0; //set counter to zero
      for(int j=j_low;j<(j_high+1);j++){
        
        id_index = j % n;
        
        // replace init_weight only if hhid at id_index has same household id
        if(hhid[id_index]==id){
          init_weight[j] = 1;
          count_size = count_size + 1;
        }
        if(count_size==size){
          break;
        }
      }
    }
  }

  return init_weight;
}

// [[Rcpp::export(sumVec)]]
IntegerVector sumVec(IntegerVector init_weight,int sizefactor){
  // loop over init_weight by a multitude of n and sizefactor
  // to create sum over households
  int n=init_weight.size()/sizefactor;
  IntegerVector init_hh(n);
  for(int i=0;i<n;i++){
    for(int j=0;j<sizefactor;j++){
      init_hh[i] = init_hh[i] + init_weight[(j*n)+i]; 
    }
  }
  
  return init_hh;
}

// [[Rcpp::export(select_equal)]]
List select_equal(IntegerVector x,int val1, int val2){
  
  int n1 = sum(x==val1);
  int n2 = sum(x==val2);
  int k1 = 0;
  int k2 = 0;
  IntegerVector out1(n1);
  IntegerVector out2(n2);
  
  for(int i=0;i<x.size();i++){
   if(x[i]==val1){
     out1[k1] = i;
     k1 = k1+1;
   }
   if(x[i]==val2){
     out2[k2] = i;
     k2 = k2+1;
   }
  }
  
  return Rcpp::List::create(Rcpp::Named("hh_1") = out1,
                     Rcpp::Named("hh_2") = out2);
  
}


// [[Rcpp::export]]
std::map<int,int> tableC(IntegerVector x){
  // Create a map
  std::map<int, int> tab;
  
  tab.clear();
  
  // Fill the map with occurrences per number.
  for (int i = 0; i < x.size(); ++i) {
    tab[ x[i] ] += 1;
  }
  return tab;
}

// [[Rcpp::export]]
IntegerVector csample_num( IntegerVector x,
                           int size,
                           bool replace,
                           NumericVector prob = NumericVector::create()){
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

// [[Rcpp::export]]
IntegerVector sample_group(IntegerVector x, IntegerVector group_x, IntegerVector group, IntegerVector group_num,bool replace){
  
  int n = group.size();
  IntegerVector out(sum(group_num));
  int out_pos =0;
  
  for(int i=0;i<n;i++){
   // sample integer values
   IntegerVector sample_i = csample_num(x[group_x==group[i]],group_num[i],replace);
   
   for(int j=0;j<sample_i.size();j++){
     out[j+out_pos] = sample_i[j];
   }
   out_pos = out_pos+sample_i.size();
  }
  return out;
}

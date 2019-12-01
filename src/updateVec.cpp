#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export(updateVecC)]]
IntegerVector updateVecC(IntegerVector &init_weight,IntegerVector &add_index, IntegerVector &remove_index, IntegerVector &hhsize, IntegerVector &hhid,int &sizefactor) {
  
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
  
  /*
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
   */
  //return add_id;
  // got through add_index and remove_index
  // and select update init_weight using hhsize and hhid
  for(int i=0;i<n_ar;i++){
    
    if(i<n_remove){
      //remove household 
      // select size and id using modulus operator
      id_index = remove_index[i] % n;
      size = hhsize[id_index];
      id = hhid[id_index];
      
      j_low = remove_index[i]-(size-1);
      j_low = std::max(remove_index[i]-id_index,j_low);
      j_high = remove_index[i]+(size-1);
      j_high = std::min(j_high,remove_index[i]+n-id_index);
      
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
      id_index = add_index[i] % n;
      size = hhsize[id_index];
      id = hhid[id_index];
      
      j_low = add_index[i]-(size-1);
      j_low = std::max(add_index[i]-id_index,j_low);
      j_high = add_index[i]+(size-1);
      j_high = std::min(j_high,add_index[i]+n-id_index);
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

// helpfunction to split vector between 0 and 1
// [[Rcpp::export]]
Rcpp::List splitVector(Rcpp::IntegerVector &x){
  
  int ones = 0;
  int zeros = 0;
  int nOnes = sum(x);
  Rcpp::IntegerVector xOnes(nOnes);
  Rcpp::IntegerVector xZeros(x.size()-nOnes);
  
  for(int i=0;i<x.size();i++){
    if(x[i]==1){
      xOnes[ones] = i;
      ones += 1;
    }else{
      xZeros[zeros] = i;
      zeros += 1;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("indexAdd") = xZeros,
                            Rcpp::Named("indexRemove") = xOnes);
  
}

// helpfunction to calcualte probabilities for each
// case
// called inside c++ function calcProbabilities()
//
// [[Rcpp::export]]
double calcCase(Rcpp::NumericVector &x){

  double probx = 0;
  if(sum(x)!=0){
    // if different from 0
    // calculate mean of x multiplied by weighted sum over sign(x)
    double meanx = std::abs(mean(x));
    //Rcpp::NumericVector signx = as<NumericVector>(sign(x));
    double wmSign = sum(as<NumericVector>(sign(x))*x)/sum(x);

    probx = meanx*wmSign;
    
  }
  
  return probx;
  
}  

// help function to calculate probabilites using x and an index matrix
// used to calculate sampling probabilities
// probabilities <- sum(x[indexMat(i,_)])
// if probabilities <=0 ==> exp(sum(x[indexMat(i,_)]))
// x = vector containing differencen
// indexMat = matrix containing indices which subset x
// initWeight = 0-1 vector
// [[Rcpp::export]]
Rcpp::List calcProbabilities(Rcpp::IntegerMatrix &indexMat, Rcpp::NumericVector &x, Rcpp::IntegerVector &indexData, Rcpp::IntegerVector &initWeight,
                             Rcpp::IntegerVector &indexAdd, Rcpp::IntegerVector &indexRemove){
  
  int nrow = indexMat.nrow();
  Rcpp::NumericVector probAdd(nrow);
  Rcpp::NumericVector helpVec(indexMat.ncol());
  Rcpp::LogicalVector negIndex(nrow);
  double sumNegatives = 0.0;
  double sumPositives = 0.0;
  
  // get probabilites for adding
  for(int i=0;i<nrow;i++){
    helpVec = x[indexMat(i,_)];
    probAdd[i] = calcCase(helpVec);
    negIndex[i] = probAdd[i]<=0;
    if(negIndex[i]){
      sumNegatives += probAdd[i];
    }else{
      sumPositives += probAdd[i];
    }
  }
  
  // get probabilities for removing
  Rcpp::NumericVector probRemove = probAdd*-1;
  
  // addjust probabilities for negative differences
  probAdd[negIndex] = exp(sumNegatives);
  probRemove[!negIndex] = exp(-1*sumPositives);
  
  // get probabilites for each index in indexData
  // considering initWeight
  int helpIndex = 0;
  int kRemove = 0;
  int kAdd = 0;
  int sizeData = indexData.size();
  int Ones = indexRemove.size();
  int Zeros = indexAdd.size();
  std::vector<double> probAdd_out(Zeros);
  std::vector<double> probRemove_out(Ones);
  std::vector<int> indexRemoveNew(Ones);
  std::vector<int> indexAddNew(Zeros);
  

  for(int i=0; i<std::max(Ones,Zeros);i++){
    helpIndex =  indexData[i % sizeData];
    if(i<Zeros && probAdd[helpIndex]>0){
      indexAddNew[kAdd] = indexAdd[i];
      probAdd_out[kAdd] = probAdd[helpIndex];
      kAdd = kAdd +1;
    }
    if(i<Ones && probRemove[helpIndex]>0){
      indexRemoveNew[kRemove] = indexRemove[i];
      probRemove_out[kRemove] = probRemove[helpIndex];
      kRemove = kRemove+1;
    }
  }
  probRemove_out.resize(kRemove);
  indexRemoveNew.resize(kRemove);
  probAdd_out.resize(kAdd);
  indexAddNew.resize(kAdd);
  
  
  return Rcpp::List::create(Rcpp::Named("probAdd") = probAdd_out,
                            Rcpp::Named("probRemove") = probRemove_out,
                            Rcpp::Named("indexRemove") = indexRemoveNew,
                            Rcpp::Named("indexAdd") = indexAddNew,
                            Rcpp::Named("nAdd") = kAdd,
                            Rcpp::Named("nRemove") = kRemove);
}



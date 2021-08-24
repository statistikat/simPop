#include <RcppArmadilloExtensions/sample.h>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <algorithm>    // std::count
#include <vector>       // std::vector
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]] 

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
    double wmSign = sum(abs(x))/sum(x);
    
    probx = meanx*wmSign;
    
  }
  
  return probx;
  
}  

// helpfunction to calcualte probabilities for each
// case
// called inside c++ function calcProbabilities()
//
// [[Rcpp::export]]
double calcCase2(Rcpp::NumericVector &x){
  
  Rcpp::NumericVector x_2 = as<NumericVector>(sign(x)) * pow(x,2.0);
  double probx = sum(x_2);
  probx = sqrt(abs(probx)) * std::copysign(1.0,probx);
  
  return(probx);
  
}  

// help function to calculate probabilites using x and an index matrix
// used to calculate sampling probabilities
// probabilities <- sum(x[indexMat(i,_)])
// if probabilities <=0 ==> exp(sum(x[indexMat(i,_)]))
// x = vector containing differencen
// indexMat = matrix containing indices which subset x
// initWeight = 0-1 vector
// [[Rcpp::export]]
Rcpp::List calcProbabilities(Rcpp::IntegerMatrix &indexMat, Rcpp::NumericVector &x, Rcpp::NumericVector &Npop, Rcpp::IntegerVector &indexData, Rcpp::IntegerVector &initWeight,
                             Rcpp::IntegerVector &indexAdd, Rcpp::IntegerVector &indexRemove, int n_add, int n_remove){
  
  int nrow = indexMat.nrow();
  Rcpp::NumericVector probAdd(nrow);
  Rcpp::NumericVector helpVec(indexMat.ncol());
  Rcpp::LogicalVector negIndex(nrow);
  double sumNegatives = 0.0;
  double sumPositives = 0.0;
  // std::cout<<"start loop";
  // get probabilites for adding
  for(int i=0;i<nrow;i++){
    
    // std::cout<<i<<"   ";
    helpVec = x[indexMat(i,_)];
    probAdd[i] = calcCase(helpVec);
    negIndex[i] = probAdd[i]<=0;
    if(negIndex[i]){
      sumNegatives += probAdd[i];
    }else{
      sumPositives += probAdd[i];
    }
  }
  // std::cout<<"finish\n";
  // get probabilities for removing
  Rcpp::NumericVector probRemove = probAdd*-1;
  
  // addjust probabilities for negative differences
  probAdd[negIndex] = exp(sumNegatives);
  probRemove[!negIndex] = exp(-1*sumPositives);
  
  // std::cout<<"start adjustment";
  if(max(x)<n_add){
    // std::cout<<"adjust\n";
    // create weighted mean between probAdd and probRemove 
    // if number of draws succeeds highest positived difference to target margins
    probAdd = (probAdd*max(x) + probRemove*(n_add-max(x)))/n_add;
  }
  if(abs(min(x))<n_remove){
    probRemove = (probRemove*abs(min(x)) + probAdd*(n_remove-abs(min(x))))/n_remove;
  }
  
  // for(int i=0;i<probAdd.size();i++){
  //   std::cout<<"probAdd"<<probAdd[i]<<"\n";
  // }
  // 
  // for(int i=0;i<probRemove.size();i++){
  //   std::cout<<"probRemove"<<probRemove[i]<<"\n";
  // }
  
  // get probabilites for each index in indexData
  // considering initWeight
  int helpIndex = 0;
  int kRemove = 0;
  int kAdd = 0;
  int sizeData = indexData.size();
  int Ones = indexRemove.size();
  int Zeros = indexAdd.size();
  
  /*
   Rcpp::NumericVector probAdd_help = probAdd[indexData][indexAdd % sizeData];
   Rcpp::NumericVector probRemove_help = probRemove[indexData][indexRemove % sizeData];
   Rcpp::LogicalVector addPos = probAdd_help>0;
   Rcpp::LogicalVector removePos = probRemove_help>0;
   Rcpp::IntegerVector indexAddNew = indexAdd[addPos];
   Rcpp::IntegerVector indexRemoveNew = indexRemove[removePos];
   Rcpp::NumericVector probAdd_out = probAdd_help[addPos];
   Rcpp::NumericVector probRemove_out = probRemove_help[removePos];
   */
  std::vector<double> probAdd_out(Zeros);
  std::vector<double> probRemove_out(Ones);
  std::vector<int> indexRemoveNew(Ones);
  std::vector<int> indexAddNew(Zeros);
  std::vector<int> count_group_add(nrow);
  std::vector<int> count_group_remove(nrow);
  
  // loop through vectors to add and remove to get number of units in each group
  for(int i=0; i<std::max(Ones,Zeros);i++){
    if(i<Zeros){ // && 
      helpIndex =  indexData[indexAdd[i] % sizeData];
      count_group_add[helpIndex]++;
    }
    if(i<Ones){ //
      helpIndex =  indexData[indexRemove[i] % sizeData];
      count_group_remove[helpIndex]++;
    }
  }
  // adjust probabilites to the population size
  for(int i=0;i<nrow;i++){
    probAdd[i] = probAdd[i]/count_group_add[i];
    probRemove[i] = probRemove[i]/count_group_remove[i];
  }
  
  // create vector with sample probabilites for whole possible sampel size
  for(int i=0; i<std::max(Ones,Zeros);i++){
    if(i<Zeros){ // && 
      helpIndex =  indexData[indexAdd[i] % sizeData];
      if(probAdd[helpIndex]>0){
        indexAddNew[kAdd] = indexAdd[i];
        probAdd_out[kAdd] = probAdd[helpIndex];
        kAdd = kAdd +1;
      }
    }
    if(i<Ones){ //
      helpIndex =  indexData[indexRemove[i] % sizeData];
      if(probRemove[helpIndex]>0){
        indexRemoveNew[kRemove] = indexRemove[i];
        probRemove_out[kRemove] = probRemove[helpIndex];
        kRemove = kRemove+1;
      }
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



//  function to update difference vector and init_weight
// difference vector is calculated as target counts - synthetic counts
// function returns new difference between target margins and current margins
// [[Rcpp::export(updateObjectiveC)]]
Rcpp::List updateObjectiveC(IntegerVector &init_weight,
                            IntegerVector &add_index,
                            IntegerVector &remove_index,
                            IntegerVector &hhsize,
                            IntegerVector &hhid,int &sizefactor,
                            IntegerMatrix &indexMat,
                            IntegerVector &indexData,
                            NumericVector diff,
                            IntegerVector &householdMargin) {
  
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
  bool house_margins_updated = false; // for internal loop
  std::unordered_set<int> indexUsed;
  int ncol = indexMat.ncol();
  IntegerVector indexVec(ncol);
  IntegerVector help_count(ncol);
  
  // std::cout<< ncol;
  // got through add_index and remove_index
  // and select update init_weight using hhsize and hhid
  // update diff accordingly
  // ~ each time an index is added or removed add or substract 1 to all marigins the index is part of
  // each household added the household margins are updated only once using house_marings_updated
  for(int i=0;i<n_ar;i++){
    
    // std::cout<<"test1\n";
    if(i<n_remove){
      // remove household if it was not already removed by other index
      if(indexUsed.find(remove_index[i])==indexUsed.end()){
        // if index is not in map start removing households which belong to index
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
            // std::cout<<"remove\n";
            // update difference vector
            // std::cout<<"test again "<<indexData[id_index]<<"\n";
            indexVec = indexMat(indexData[id_index],_);
            for(int d=0; d<ncol; d++){
              // std::cout<<"test d"<<indexVec[d]<<"\n";
              // std::cout<<"test again"<<diff[indexVec[d]]<<"\n";
              if(householdMargin[indexVec[d]]==0){
                // difference vector is target - synthetic
                // removing an index from synthetic -> diff is added by 1
                diff[indexVec[d]]++; 
              }else if(householdMargin[indexVec[d]]==1 && house_margins_updated==false){
                diff[indexVec[d]]++;
              }
            }
            house_margins_updated = true; // after the first person in a household was removed dont update household tables again
            indexUsed.insert(j); // save index in map so it is not used again
          }
          if(count_size==size){
            break;
          }
        }
        // set  house_margins_updated to false again for next round
        house_margins_updated = false;
      }
    }
    
    if(i<n_add){
      // add household
      // add household if it was not already added by other index
      if(indexUsed.find(add_index[i])==indexUsed.end()){
        // if index is not in map start adding households which belong to index
        // select size and id using modulus operator
        // std::cout<<"add\n";
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
            
            // update difference vector
            indexVec = indexMat(indexData[id_index],_);
            for(int d=0; d<ncol; d++){
              if(householdMargin[indexVec[d]]==0){
                // difference vector is target - synthetic
                // adding an index from synthetic -> diff is reduced by 1
                diff[indexVec[d]]--;
              }else if(householdMargin[indexVec[d]]==1 && house_margins_updated==false){
                diff[indexVec[d]]--;
              }
            }
            house_margins_updated = true; // after the first person in a household was removed dont update household tables again
            indexUsed.insert(j); // save index in map so it is not used again
          }
          if(count_size==size){
            break;
          }
        }
        // set  house_margins_updated to false again for next round
        house_margins_updated = false;
      }
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("init_weight") = init_weight,
                            Rcpp::Named("diff_new") = diff);
}

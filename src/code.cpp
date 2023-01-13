#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export(rng = false)]]
int find_max_from_sel_fast(NumericVector x, IntegerVector which){
  int val = 1;
  double max = x[which[0]-1];
  int which_length = which.length();
  for(int i = 1; i < which_length; i++){
    int w = which[i];
    if(IntegerVector::is_na(w)){;
      break;
    }else{
      if(x[w-1] > max){
        max = x[w-1];
        val = i + 1;
      }
    }
  }
  return val;
}


//[[Rcpp::export(rng = false)]]
int find_min_from_sel_fast(NumericVector x, IntegerVector which){
  int val = 1;
  double min = x[which[0]-1];
  int which_length = which.length();
  for(int i = 1; i < which_length; i++){
    int w = which[i];
    if(IntegerVector::is_na(w)){;
      break;
    }else{
      if(x[w-1] < min){
        min = x[w-1];
        val = i + 1;
      }
    }
  }
  return val;
}

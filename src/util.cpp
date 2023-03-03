#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export(rng = FALSE)]]
NumericVector cumz(NumericVector x){
  int size = x.size();
  NumericVector res(size);
  if(size == 0){
    return res;
  }else if(size == 1){
    res[0] = NA_REAL;
    return res;
  }
  double m = x[0];
  double msq = 0;
  res[0] = NA_REAL;

  for(int i = 2; i <= size; i++){
    double delta = x[i-1] - m;
    m += delta / i;
    double delta2 = x[i-1] - m;
    msq = (msq * (i - 1) + delta * delta2) / i;
    res[i-1] = m / sqrt(msq / (i-1));
  }
  return res;
}


// [[Rcpp::export(rng = FALSE)]]
List cumz_which_abs_max(NumericVector x){
  int size = x.size();
  int max_idx = 0;
  double max = -std::numeric_limits<float>::infinity();
  if(size == 0){
    return max_idx;
  }else if(size == 1){
    return 1;
  }
  double m = x[0];
  double msq = 0;

  for(int i = 2; i <= size; i++){
    double delta = x[i-1] - m;
    m += delta / i;
    double delta2 = x[i-1] - m;
    msq = (msq * (i - 1) + delta * delta2) / i;
    double val = abs(m / sqrt(msq / (i-1)));
    if(val > max){
      max = val;
      max_idx = i;
    }
  }
  return List::create(Named("max") = max , Named("idx") = max_idx);
}



/*** R
timesTwo(42)

cum_z_stat2 <- function(x){
  out <- rep(NA_real_, length(x))
  m <- x[1]
  msq <- 0
  out[1] <- NA
  for(idx in seq_along(x)[-1]){
    delta <- x[idx] - m
    m <- m + delta / idx
    delta2 <- x[idx] - m
    msq <- (msq * (idx-1) + delta * delta2) / idx

    out[idx] <- m / sqrt(msq / (idx-1))
  }
  out
}


x <- rnorm(16000)
bench::mark(
  which.max(abs(cum_z_stat2(x))),
  which.max(abs(cumz(x))),
  cumz2(x),
  check = TRUE
)

*/

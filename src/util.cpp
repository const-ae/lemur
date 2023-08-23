#include <RcppArmadillo.h>
using namespace Rcpp;


// // [[Rcpp::export(rng = FALSE)]]
// NumericVector cumz(NumericVector x){
//   int size = x.size();
//   NumericVector res(size);
//   if(size == 0){
//     return res;
//   }else if(size == 1){
//     res[0] = NA_REAL;
//     return res;
//   }
//   double m = x[0];
//   double msq = 0;
//   res[0] = NA_REAL;
//
//   for(int i = 2; i <= size; i++){
//     double delta = x[i-1] - m;
//     m += delta / i;
//     double delta2 = x[i-1] - m;
//     msq = (msq * (i - 1) + delta * delta2) / i;
//     res[i-1] = m / sqrt(msq / (i-1));
//   }
//   return res;
// }


// [[Rcpp::export(rng = FALSE)]]
List cumz_which_abs_max(NumericVector x, int min_neighborhood_size){
  int size = x.size();
  min_neighborhood_size = std::min(min_neighborhood_size, size);
  int max_idx = 0;
  double max = -std::numeric_limits<float>::infinity();
  int sign = +1;
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
    if(i >= min_neighborhood_size - 1 && val > max){
      max = val;
      sign = m < 0 ? -1 : +1;
      max_idx = i;
    }
  }
  return List::create(Named("max") = sign * max , Named("idx") = max_idx);
}


// [[Rcpp::export(rng = FALSE)]]
List cum_brls_which_abs_max(const NumericVector y, const arma::mat& X, const IntegerVector group,
                            const arma::rowvec& contrast, const double penalty, int min_neighborhood_size){
  int size = y.size();
  min_neighborhood_size = std::min(min_neighborhood_size, size);

  int max_idx = 0;
  double max_val = -std::numeric_limits<float>::infinity();
  int sign = +1;
  if(size == 0){
    return  List::create(Named("index") = max_idx, Named("max") = max_val);
  }else if(size == 1){
    return  List::create(Named("index") = 1, Named("max") = y);
  }
  int k = X.n_cols;
  int g = max(group);
  arma::colvec m(g);
  IntegerVector count(g);
  arma::mat X_act(g, k);
  arma::mat gamma = 1/penalty * arma::eye(k, k);
  arma::colvec beta(k);
  int n_obs = 0;
  double se_pre = 0;

  for(int i = 0; i < size; ++i){
    double yi = y[i];
    arma::colvec xi = X.row(i).t();
    int gi = group[i] - 1;
    double delta_m = 1/(count[gi] + 1.0) * yi - (1 - count[gi] / (count[gi] + 1.0)) * m(gi);
    m(gi) += delta_m;
    count(gi) += 1;
    if(count(gi) == 1){
      X_act.row(gi) = xi.t();
      n_obs += 1;
      gamma -= (gamma * xi * xi.t() * gamma) / arma::as_scalar(1 + xi.t() * gamma * xi);
      se_pre = arma::as_scalar(contrast * gamma * contrast.t());
      beta += gamma * xi * (m[gi] - xi.t() * beta);
    }else{
      beta += gamma * (xi * delta_m);
    }
    if(n_obs > k){
      double rss = std::max(1e-6, arma::accu(arma::pow(m - X_act * beta, 2)));
      double se = sqrt(se_pre * rss / (n_obs - k));
      double t_stat = arma::as_scalar(contrast * beta) / se;
      if(i >= min_neighborhood_size-1 && abs(t_stat) > max_val){
        max_val = abs(t_stat);
        sign = t_stat < 0 ? -1 : +1;
        max_idx = i + 1;
      }
    }
  }
  return List::create(Named("idx") = max_idx, Named("max") = sign * max_val);
}


// [[Rcpp::export(rng = FALSE)]]
IntegerVector count_neighbors_fast(NumericMatrix knn_mat, IntegerVector indices){
  int n_cells = knn_mat.nrow();
  int knn = knn_mat.ncol();
  int n = indices.length();
  IntegerVector counter(n_cells);
  for(int i = 0; i < n; ++i){
    for(int k = 0; k < knn; ++k){
      counter(knn_mat(indices(i)-1, k)-1) += 1;
    }
  }
  return counter;
}


/*** R

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

cum_z_stat2(x[1:10])
idx <- 10
.x <- x[seq_len(idx)]
msq <- sum((mean(.x) - .x)^2)
mean(.x) / sqrt(msq / idx / (idx-1))
mean(.x) / (sd(.x) / sqrt(idx))

x <- rnorm(16000)
cumz_which_abs_max(x)
max(abs(cum_z_stat2(x)), na.rm = TRUE)
bench::mark(
  which.max(abs(cum_z_stat2(x))),
  which.max(abs(cumz(x))),
  cumz2(x),
  check = TRUE
)

*/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cstdint>
#include <cstring>
#include "harmony_defines.h"

// [[Rcpp::export]]
void set_harmony_Zcorr(SEXP harmonyObj, const arma::Mat<float>& Z_new) {
    Rcpp::XPtr<harmony> ptr(harmonyObj);
    ptr->Z_corr = Z_new;
}

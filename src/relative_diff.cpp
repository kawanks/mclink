// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Convert sparse column in a dense vector with: -1, 0 and 1
arma::Col<int> convert_sign_vector(const arma::sp_vec &v, arma::uword n_rows) {
  arma::Col<int> result(n_rows, arma::fill::zeros);
  
  for (arma::uword i = 0; i < v.n_nonzero; ++i) {
    arma::uword row = v.row_indices[i];
    double val = v.values[i];
    result(row) = (val < 0) ? -1 : (val > 0 ? 1 : 0);
  }
  
  return result;
}

// Compare two vectors and calculate the proportion of equalities ignoring -1
double hamming_proportion(const arma::Col<int> &a, const arma::Col<int> &b) {
  arma::uword n = a.n_elem;
  arma::uword matches = 0, valid = 0;
  
  for (arma::uword i = 0; i < n; ++i) {
    if (a[i] < 0 || b[i] < 0) continue;
    valid++;
    if (a[i] == b[i]) matches++;
  }
  
  return (valid > 0) ? static_cast<double>(matches) / valid : -1.0;
}

// [[Rcpp::export]]
arma::mat relative_diff(const arma::sp_mat &M) {
  arma::uword ncols = M.n_cols;
  arma::uword nrows = M.n_rows;
  arma::mat result(ncols, ncols);
  
  // Convert all sparse columns in dense columns 
  std::vector<arma::Col<int>> sign_cols(ncols);
  for (arma::uword i = 0; i < ncols; ++i) {
    sign_cols[i] = convert_sign_vector(M.col(i), nrows);
  }
  
  for (arma::uword i = 0; i < ncols; ++i) {
    for (arma::uword j = i + 1; j < ncols; ++j) {
      double dist = 1 - hamming_proportion(sign_cols[i], sign_cols[j]);
      result(i, j) = dist;
      result(j, i) = dist;
    }
  }
 
  return result;
}


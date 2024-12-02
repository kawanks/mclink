// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::sp_mat relative_diff(arma::sp_mat &M_transpose)
{
  // Transpose the genomic matrix for operating in rows
  arma::sp_mat M = M_transpose.t();
  
  size_t rows = M.n_rows;
  size_t cols = M.n_cols;
  arma::sp_mat relative_diff_matrix(rows, rows);
  
  for(arma::sp_mat::const_iterator it = M.begin(); it != M.end(); ++it)
  {
    size_t row = it.row();
    size_t col = it.col();
    
    size_t cont = 0;
    
    for(size_t k = col + 1; k < cols; ++k)
    {
      if (M(col + 1) <= 0) {
        cont += 1;
      }
    }
    
    
  }
  
  /*
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = i+1; j < rows; ++j) {
      double count = 0.0;
      double denom = 0.0;
      
      for (size_t k = 0; k < cols; ++k) {
        if ((M(i, k) == 0 && M(j, k) == 0) || (M(i, k) > 0 && M(j, k) > 0)) {
          count += 1;
        }
        if (M(i, k) >= 0 && M(j, k) >= 0) {
          denom += 1;
        }
      }
      
      if (denom > 0) {
        relative_diff_matrix(i, j) = relative_diff_matrix(j, i) = 1.0 - (count / denom); // Tells how divergent they are
      }
    }
  }*/
  
  return relative_diff_matrix;
}

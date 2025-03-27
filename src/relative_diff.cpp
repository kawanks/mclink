// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::sp_mat relative_diff(arma::sp_mat &M)
{
  // Get the dimensions of the matrix
  size_t cols = M.n_cols;
  size_t rows = M.n_rows;
  
  /*
   * The matrix of relative differences is defined as sparse matrix just because
   * the operations between a sparse matrix type and a dense matrix type cost 
   * more, in terms of complexity, than the operations between two matrices of
   * sparse type, even that one of them is a dense matrix.
   */
  arma::sp_mat relative_diff_matrix(cols, cols);
  
  for(size_t i = 0; i < cols; i++)
  {
    for(size_t k = i + 1; k < cols; k++)
    {
      size_t diff = 0;
      
      for(arma::sp_mat::const_col_iterator it = M.begin_col(i); it != M.end_col(i); ++it)
      {
        size_t row = it.row();
        
        if(M(row, k) == 0) diff++;
        else if(M(row, k) < 0) rows--;
        
      }
      
      for(arma::sp_mat::const_col_iterator it = M.begin_col(k); it != M.end_col(k); ++it)
      {
        size_t row = it.row();
        
        if(M(row, k) == 0) diff++;
        else if(M(row, k) < 0) rows--;
      }
      
      relative_diff_matrix(i, k) = (double) diff/rows;
    }
    
  }
  
  return relative_diff_matrix;
}

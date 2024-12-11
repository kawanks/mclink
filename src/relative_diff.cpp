// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
arma::sp_mat relative_diff(arma::sp_mat &M)
{
  // Get the dimensions of the matrix
  size_t rows = M.n_rows;
  size_t cols = M.n_cols;
  
  /*
   * The matrix of relative differences is defined as sparse matrix just because
   * the operations between a sparse matrix type and a dense matrix type cost 
   * more, in terms of complexity, than the operations between two matrices of
   * sparse type, even that one of them is a dense matrix.
   */
  arma::sp_mat relative_diff_matrix(rows, rows);
  
  Progress p(rows, true);
  for(size_t i = 0; i < rows - 1; i++)
  {
    if (Progress::check_abort() )
      return NULL;
    
    for(size_t k = i + 1; k < rows; k++)
    {
      size_t diff = 0;
      
      for(arma::sp_mat::const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it)
      {
        size_t col = it.col();
        
        if(M(k, col) <= 0) diff++;
      }
      
      for(arma::sp_mat::const_row_iterator it = M.begin_row(k); it != M.end_row(k); ++it)
      {
        size_t col = it.col();
        
        if(M(i, col) <= 0) diff++;
      }
      
      relative_diff_matrix(i, k) = relative_diff_matrix(k, i) = (double)diff/cols;
    }
    
    p.increment();
    
    Rprintf("Linha %zu\n", i);
  }
  
  return relative_diff_matrix;
}

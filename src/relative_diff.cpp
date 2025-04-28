// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::sp_mat relative_diff(arma::sp_mat &M)
{
  // Get the dimensions of the matrix
  size_t cols = M.n_cols, diff, rows;
  
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
      diff = 0;
      rows = M.n_rows;
      
      for(arma::sp_mat::const_col_iterator it = M.begin_col(i); it != M.end_col(i); ++it)
      {
        size_t row = it.row();
        
        if(M(row, i) < 0 || M(row, k) < 0) rows--;
        else if(M(row, k) == 0) diff++;
      }
      
      for(arma::sp_mat::const_col_iterator it = M.begin_col(k); it != M.end_col(k); ++it)
      {
        size_t row = it.row();
        
        if(M(row, i) == 0 && M(row, k) < 0) rows--;
        else if(M(row, i) == 0 && M(row, k) > 0) diff++;
      }
      
      relative_diff_matrix(i, k) = (double) diff/rows;
    }
    
  }
  
  return relative_diff_matrix;
}


// [[Rcpp::export]]
arma::sp_mat shuffle_relative_diff(arma::sp_mat &M)
{
  arma::mat shuffled = arma::shuffle(M.as_dense(), 1);
  arma::sp_mat sp_shuffled = arma::sp_mat(shuffled);
  
  return relative_diff(sp_shuffled);
}


//[[Rcpp:export]]
double cal_var(arma::sp_mat &M)
{
  return arma::var(arma::nonzeros(M));
}


// [[Rcpp::export]]
arma::vec simulate_panmixia(arma::sp_mat &dados, size_t iterations)
{
  arma::sp_mat e_i;
  
  for(size_t i = 0 ; i < iterations; ++i){
    e_i = shuffle_relative_diff(dados);
  }
  
  return arma::nonzeros(e_i);
}







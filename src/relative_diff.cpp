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
arma::sp_mat shuffle_relative_diff(arma::sp_mat M)
{
  // Shuffle the matrix M passed as parameter
  size_t n_cols = M.n_cols, n_rows = M.n_rows;
  arma::vec values(n_cols);
  arma::sp_mat shuffled(n_rows, n_cols);
  
  for(size_t i = 0; i < n_rows; ++i){
    values.zeros();
    for(arma::sp_mat::const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it){
      values[it.col()] = *it;
    }
    
    values = arma::shuffle(values); 
    
    for(size_t k = 0; k < n_cols; ++k){
      if(values[k] != 0){
        shuffled(i, k) = values[k];
      }
    }
  }
  
  // return the relative distances of the matrix shuffled.
  return relative_diff(shuffled);
}


// [[Rcpp::export]]
arma::vec simulate_panmixia(arma::sp_mat &dados, size_t iterations)
{
  // Simulate the panmixia shuffling the data matrix a number 'iterations'
  // times.
  arma::vec means(iterations), var(iterations), parameters(2);
  arma::sp_mat e_i(arma::size(dados));
  size_t sample_size;
  
  for(size_t i = 0 ; i < iterations; ++i){
    e_i = shuffle_relative_diff(dados);
    
    means(i, 0) = arma::mean(arma::nonzeros(e_i)); // The mean of the sample
    
    var[i] = arma::var(arma::nonzeros(e_i)); // The variance of the sample
    
    std::cout<< i*100/iterations << "%\n";
  }
  
  std::cout<< "100%\n";
  
  sample_size = arma::vec(arma::nonzeros(e_i)).n_elem;
  
  parameters(0) = arma::sum(means) / iterations; // mean of all elements
  parameters(1) = arma::sum(var) * (sample_size - 1) / (iterations * sample_size - 1); // variance of all elements
  
  // return the mean and the variance of all simulations.
  return parameters;
}
  

// [[Rcpp::export]]
arma::mat simulate_panmixia_vec(arma::sp_mat &dados, size_t iterations)
{
  
  arma::mat elements((dados.n_cols*dados.n_cols - (dados.n_cols))/2, iterations);
  arma::sp_mat e_i(arma::size(dados));
  
  for(size_t i = 0 ; i < iterations; ++i){
    e_i = shuffle_relative_diff(dados);
    
    elements.col(i) = arma::vec(arma::nonzeros(e_i)); // The mean of the sample
    
    std::cout<< i*100/iterations << "%\n";
  }
  
  // return all elements from all the simulations.
  return elements;
}
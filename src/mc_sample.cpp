// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>
#include <string>
#include <algorithm>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "relative_diff.h"


// [[Rcpp::export(name = ".mc_sample_rows")]]
std::vector<int> mc_sample_rows(const std::vector<int>& chr, 
                                const std::vector<int>& pos,
                                const size_t sample_size = 1000,
                                const int min_distance = 1000){
  
  int size = pos.size();
  std::vector<int> rows(size);
  std::iota(rows.begin(), rows.end(), 0);
  
  std::vector<int> sample;
  sample.reserve(sample_size);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  
  while(!rows.empty() && sample.size() < sample_size){
    std::uniform_int_distribution<int> unif(0, rows.size() - 1);
    int rpos = unif(gen);
    int sampled_row = rows[rpos];
    
    sample.push_back(sampled_row);
    
    int center = pos[sampled_row];

    int sampled_chr = chr[sampled_row];
    
    auto lower = std::lower_bound(rows.begin(), rows.end(), center - min_distance, 
                                  [&](int row, int value){
                                    if(chr[row] != sampled_chr) return chr[row] < sampled_chr;
                                    return pos[row] < value;
                                  });
    
    auto upper = std::upper_bound(rows.begin(), rows.end(), center + min_distance,
                                  [&](int value, int row){
                                    if(chr[row] != sampled_chr) return chr[row] > sampled_chr;
                                    return pos[row] > value;
                                  });
    
    rows.erase(lower, upper);
  }


  // Since R starts the index of the vectors in 1 is needed to add 1 on the sample vector to return
  std::for_each(sample.begin(), sample.end(), [](int &n){
          n += 1;
          });

  return sample;
}

arma::sp_mat mc_sample_matrix(const Rcpp::List data, const size_t sample_size, const int min_distance) {
  arma::sp_mat M = Rcpp::as<arma::sp_mat>(data[0]);
  std::vector<int> chrs = data[1];
  std::vector<int> pos = data[2];

  std::vector<int> rows = mc_sample_rows(chrs, pos, sample_size, min_distance);
  arma::uword n_sample = rows.size();
  arma::uword n_cols = M.n_cols;
  
  std::vector<arma::uword> row_inds;
  std::vector<arma::uword> col_inds;
  std::vector<double> values;
  
  // Optimization: allocate a little more memory than necessary.
  row_inds.reserve(M.n_nonzero);  
  col_inds.reserve(M.n_nonzero);
  values.reserve(M.n_nonzero);
  
  for (arma::uword i = 0; i < n_sample; ++i) {
    arma::uword row = rows[i];
    for (arma::sp_mat::const_row_iterator it = M.begin_row(row); it != M.end_row(row); ++it) {
      row_inds.push_back(i);          
      col_inds.push_back(it.col());   
      values.push_back(*it);          
    }
  }
  
  arma::umat locations(2, values.size());
  
  std::copy(row_inds.begin(), row_inds.end(), locations.row(0).begin());
  std::copy(col_inds.begin(), col_inds.end(), locations.row(1).begin());
  
  // Construct the sparse matrix directly
  arma::sp_mat sample_matrix(
      locations,
      arma::vec(values),
      n_sample, n_cols,
      true, false
  );
  
  return sample_matrix;
}


arma::sp_mat mc_shuffle_matrix(const arma::sp_mat &M) {
  const arma::uword n_rows = M.n_rows;
  const arma::uword n_cols = M.n_cols;
  
  std::vector<arma::uword> row_inds;
  std::vector<arma::uword> col_inds;
  std::vector<double> values;
  
  row_inds.reserve(M.n_nonzero);
  col_inds.reserve(M.n_nonzero);
  values.reserve(M.n_nonzero);
   
  std::random_device rd;
  std::mt19937 gen(rd());
   
  for (arma::uword i = 0; i < n_rows; ++i) {
    std::vector<arma::uword> original_cols;
    std::vector<double> original_vals;
     
    for (arma::sp_mat::const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it) {
      original_cols.push_back(it.col());
      original_vals.push_back(*it);
    }
     
    arma::uword nnz = original_cols.size();
    if (nnz == 0) continue;
     
    // Generate shuffled columns 
    std::vector<arma::uword> shuffled_cols(n_cols);
    std::iota(shuffled_cols.begin(), shuffled_cols.end(), 0);
    std::shuffle(shuffled_cols.begin(), shuffled_cols.end(), gen);
    shuffled_cols.resize(nnz);
     
    // Save new triple (line, shuffled column and original value)
    for (arma::uword j = 0; j < nnz; ++j) {
      row_inds.push_back(i);
      col_inds.push_back(shuffled_cols[j]);
      values.push_back(original_vals[j]);
    } 
  } 
  
  arma::umat locations(2, values.size());
  std::copy(row_inds.begin(), row_inds.end(), locations.row(0).begin());
  std::copy(col_inds.begin(), col_inds.end(), locations.row(1).begin());
   
  arma::vec val_vec(values);
   
  return arma::sp_mat(locations, val_vec, n_rows, n_cols);
} 


/*
 *  Esta função retorna uma matriz com duas colunas, a primeira coluna são
 *  as distâncias relativas observadas e a segunda coluna são as distâncias
 *  relativas de uma população panmítica.
*/
// [[Rcpp::export]]
Rcpp::List simulate_panmixia(Rcpp::List data, const int iterations = 1000, const size_t sample_size = 1000, const int min_distance = 1000){
  Rcpp::S4 geno = data[0];
  int n_cols = Rcpp::IntegerVector(geno.slot("Dim"))[1];

  std::vector<double> original_all;
  std::vector<double> shuffled_all;
  std::vector<double> variances(iterations);
  Progress p(iterations, TRUE);
  
  original_all.reserve(iterations * n_cols);
  shuffled_all.reserve(iterations * n_cols);
  
  for (int i = 0; i < iterations; ++i) {
    p.increment();
    arma::sp_mat sample_matrix = mc_sample_matrix(data, sample_size, min_distance);
    arma::sp_mat shuffled_matrix = mc_shuffle_matrix(sample_matrix);
    
    arma::vec original = relative_diff(sample_matrix);
    arma::vec shuffled = relative_diff(shuffled_matrix);
    variances[i] = arma::var(original) / arma::var(shuffled);
    
    original_all.insert(original_all.end(), original.begin(), original.end());
    shuffled_all.insert(shuffled_all.end(), shuffled.begin(), shuffled.end());
  }
 
  Rcpp::List result = Rcpp::List::create(original_all, shuffled_all, variances);
  
  return result;
}

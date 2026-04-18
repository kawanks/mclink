// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <vector>
#include <bitset>


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
arma::vec relative_diff(const arma::sp_mat &M) {
  arma::uword ncols = M.n_cols;
  arma::uword nrows = M.n_rows;
  std::vector<double> result_vec;
  result_vec.reserve(ncols * (ncols - 1) / 2); // Reserve 
  
  // Convert all sparse columns in dense columns 
  std::vector<arma::Col<int>> sign_cols(ncols);
  for (arma::uword i = 0; i < ncols; ++i) {
    sign_cols[i] = convert_sign_vector(M.col(i), nrows);
  }
  
  for (arma::uword i = 0; i < ncols; ++i) {
    for (arma::uword j = i + 1; j < ncols; ++j) {
      double dist = 1 - hamming_proportion(sign_cols[i], sign_cols[j]);
      result_vec.push_back(dist);
    }
  }
 
  return arma::vec(result_vec);
}


std::vector<double> sparse_dot_product(arma::sp_mat::const_row_iterator it_A,
                          arma::sp_mat::const_row_iterator end_A,
                          arma::sp_mat::const_row_iterator it_B,
                          arma::sp_mat::const_row_iterator end_B){
  std::vector<double> vec(2);
  double k_ij = 0.0;
  double missing = 0.0;
  
  while (it_A != end_A && it_B != end_B) {
    if (it_A.col() == it_B.col()) {
      
      if((*it_A) == -1 || (*it_B) == -1){
        ++it_A;
        ++it_B;
         ++missing;
         continue;
      }
      
      ++it_A;
      ++it_B;
      k_ij++;
    } else if (it_A.col() < it_B.col()) {
      if((*it_A) == -1) ++missing;
      ++it_A;
    } else {
      if((*it_B) == -1) ++missing;
      ++it_B;
    }
  }
  
  vec[0] = k_ij;
  vec[1] = missing;
  
  return vec;
}


Rcpp::List spmat2bitset(const arma::sp_mat &M){
   int snps = M.n_rows; 
   int cols = M.n_cols; 
   int parts = (cols + 63) / 64;
   
   std::vector<uint64_t> G(parts * snps, 0);
   std::vector<uint64_t> Mk(parts * snps, 0);
   
   for(int i = 0; i < snps; i++){
     for(int j = 0; j < cols; j++){
       int part = j / 64;
       int bit  = j % 64;
       Mk[i*parts + part] |= (1ULL << bit);
     }
   }
   
   for(int i = 0; i < snps; ++i){
     
     for(auto it = M.begin_row(i); it != M.end_row(i); ++it){
       
       int j = it.col();
       double val = *it;
       
       int part = j / 64;
       int bit  = j % 64;
       
       if(val == 1){
         G[i*parts + part] |= (1ULL << bit);
       }
       
       if(val == -1){
         // remover da máscara
         Mk[i*parts + part] &= ~(1ULL << bit);
       }
     }
   }
   
   return Rcpp::List::create(
     Rcpp::_["G"] = G,
     Rcpp::_["M"] = Mk,
     Rcpp::_["snps"] = snps,
     Rcpp::_["parts"] = parts
   );
}

// [[Rcpp::export]]
Rcpp::List ld_decay(Rcpp::List data, 
                    double maf = 0.05,
                    int bin_size = 1000,
                    int max_dist = 100000){
  
  arma::sp_mat M = Rcpp::as<arma::sp_mat>(data[0]); 
  std::vector<int> pos = Rcpp::as<std::vector<int>>(data[2]);
  
  Rcpp::List Bitset = spmat2bitset(M);
  std::vector<uint64_t> G = Rcpp::as<std::vector<uint64_t>>(Bitset[0]);
  std::vector<uint64_t> Mk = Rcpp::as<std::vector<uint64_t>>(Bitset[1]); 
  
  int snps = Rcpp::as<int>(Bitset[2]);
  int parts = Rcpp::as<int>(Bitset[3]);
  
  int n_bins = max_dist / bin_size;
  std::vector<double> sum_r2(n_bins, 0.0);
  std::vector<int> count(n_bins, 0);
  
  for (int i = 0; i < snps; i++) {
    for (int j = i + 1; j < snps; j++) {
      int dist = pos[j] - pos[i];
      if (dist > max_dist) break;
      
      int kij = 0;
      int ki  = 0;
      int kj  = 0;
      int nij = 0;
    
      for (int w = 0; w < parts; w++) {
      
        uint64_t gi = G[i * parts + w];
        uint64_t gj = G[j * parts + w];
      
        uint64_t mi = Mk[i * parts + w];
        uint64_t mj = Mk[j * parts + w];
      
        uint64_t valid = mi & mj;
        if (valid == 0) continue;
      
        // Verify if is possible use the popcount standart function of C++
        nij += __builtin_popcountll(valid);
      
        kij += __builtin_popcountll(gi & gj & valid);
        ki  += __builtin_popcountll(gi & valid);
        kj  += __builtin_popcountll(gj & valid);
      }
    
      if (nij < 10) continue;
      
      double pi = (double) ki / nij;
      double pj = (double) kj / nij;
    
      // Filter maf
      if (pi < maf || pi > (1.0 - maf)) continue;
      if (pj < maf || pj > (1.0 - maf)) continue;
    
      double pij = (double) kij / nij;
  
      double D = pij - pi * pj;
      double den = pi * (1.0 - pi) * pj * (1.0 - pj);
      
      if (den <= 0) continue;
      
      double r2 = (D * D) / den;
      
      int bin = dist / bin_size;
      if (bin < n_bins) {
        sum_r2[bin] += r2;
        count[bin] += 1;
      }
    }
  }
 
  std::vector<double> mean_r2(n_bins, NA_REAL);
  std::vector<int> distance(n_bins);
  
  for (int b = 0; b < n_bins; b++) {
    distance[b] = b * bin_size;
    
    if (count[b] > 0) {
      mean_r2[b] = sum_r2[b] / count[b];
    }
  }
  
  return Rcpp::List::create(
    Rcpp::_["distance"] = distance,
    Rcpp::_["mean_r2"] = mean_r2,
    Rcpp::_["count"] = count
  );
}

/*
std::vector<double> ld_r2(Rcpp::List data, double maf = 0.05){
  arma::sp_mat M = Rcpp::as<arma::sp_mat>(data[0]); 
  std::vector<int> pos = data[2];
  
  Rcpp::List Bitset = spmat2bitset(M);
  std::vector<uint64_t> G = Bitset[0];
  std::vector<uint64_t> Mk = Bitset[1];
  
  int snps = Rcpp::as<int>(Bitset[2]);
  int parts = Rcpp::as<int>(Bitset[3]);
  
  int count = 0;
  
  for(int i = 0; i < snps; i++){
    for(int j = 0; j < parts; j++){
      uint64_t valid = G[i*snps + j] & Mk[i*snps + j];
      count += __builtin_popcountll(valid);
    }
  }
  
  arma::vec k = M * arma::ones<arma::vec>(M.n_cols);
  arma::uword rows = M.n_rows;
  arma::uword cols = M.n_cols;
  
  std::vector<double> r2_values;
  std::vector<double> k_ij(2);
  
  for(size_t i = 0; i < rows; ++i){
    double pi = k[i] / cols;
    if(pi < maf || (1 - pi) < maf) continue;
    
    for(size_t j = i + 1; j < rows; ++j){
      double pj = k[j] / cols;
      if(pj < maf || (1 - pj) < maf) continue;
      
      k_ij = sparse_dot_product(M.begin_row(i), M.end_row(i), M.begin_row(j), M.end_row(j));
      double pij = k_ij[0] / (cols - k_ij[1]);
      
      double den = (pi*(1-pi)*pj*(1-pj));
      if(den <= 0) continue;
      
      double r2 = ((pij - pi*pj)*(pij - pi*pj)) / den;
      
      r2_values.push_back(r2);
    }
  }
  std::vector<double> r2_values;
  return r2_values;
}
 */


























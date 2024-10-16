#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix relative_diff(Rcpp::NumericMatrix& M)
{
  int rows = M.rows();
  Rcpp::NumericMatrix relative_diff_matrix(rows, rows);
  
  for (int i = 0; i < rows; ++i) {
    for (int j = i+1; j < rows; ++j) {
      double count = 0.0;
      double denom = 0.0;
      
      for (size_t k = 0; k < M.row(i).size(); ++k) {
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
  }
  
  return relative_diff_matrix;
}

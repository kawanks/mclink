#ifndef UTILS_H
#define UTILS_H

#include <armadillo>

arma::vec relative_diff(const arma::sp_mat &M);
std::vector<double> ld_r2(Rcpp::List data, double maf = 0.05);

#endif
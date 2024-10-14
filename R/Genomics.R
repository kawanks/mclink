library(Rcpp)

sourceCpp('../src/relative_diff.cpp')

dados <- read.csv('../data/dados.csv', header=F)
matrix <- as.matrix(dados)
matrix

relative_diff(matrix)

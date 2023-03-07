// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; 
using namespace arma;

//' sum from negative binomial density function
//'
//' Probability density function of the negative binomial distribution (written in C++)
//'
//' @param mat dgCMatrix expression counts
//' @param s numeric scaling factor
//' @param x numeric expression for reference profile
//' @param bg numeric background level
//' @param size_dnb int Dispersion parameter
//'
//' @return rowSums for matrix of densities
//' @useDynLib InSituType, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @exportPattern "^[[:alpha:]]+" 
//' @export
//' @examples
//' data("ioprofiles")
//' data("mini_nsclc")
//' bg <- Matrix::rowMeans(mini_nsclc$neg)
//' genes <- intersect(dimnames(mini_nsclc$counts)[[2]], dimnames(ioprofiles)[[1]])
//' mat <- mini_nsclc$counts[, genes]
//' x <- ioprofiles[genes, 1]
//' bgsub <- pmax(sweep(mat, 1, bg, "-"), 0)
//' s <- Matrix::rowSums(bgsub) / sum(x)
//' s[s <= 0] <- Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum(x)
//' lls(as(mat, "dgCMatrix"), s, x, bg, 10)
// [[Rcpp::export]]
Rcpp::NumericVector
lls(arma::sp_mat& mat, arma::vec& s, arma::vec& x, arma::vec& bg, int& size_dnb) {
  Rcpp::NumericVector res(s.n_rows);
  arma::vec::const_iterator x_iter = x.begin();
  for(; x_iter != x.end(); ++x_iter) {
    arma::vec::const_iterator s_iter = s.begin();
    arma::vec::const_iterator bg_iter = bg.begin();
    for(; s_iter != s.end(); ++s_iter) {
      double yhat = (*s_iter) * (*x_iter) + (*bg_iter);
      int i = s_iter - s.begin();
      int j = x_iter - x.begin();
      res(i) += R::dnbinom_mu(mat(i, j), size_dnb, yhat, 1);
      ++bg_iter;
    }
  }
  return res;
}


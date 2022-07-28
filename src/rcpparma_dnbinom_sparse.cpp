// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; 
using namespace arma;

//' negative binomial density function
//'
//' Probability density function of the negative binomial distribution (written in C++)
//'
//' @param x dgCMatrix of quantiles
//' @param mu dgCMatrix containing Mean of the distribution 
//' @param size_dnb Dispersion parameter
//'
//' @return dgCMatrix of densities
//' @useDynLib MLEcell, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @exportPattern "^[[:alpha:]]+" 
//' @export
// [[Rcpp::export]]
arma::sp_mat dnbinom_sparse(arma::sp_mat& x, arma::sp_mat& mu, int& size_dnb) {
  arma::sp_mat res = arma::sp_mat(size(x));
  arma::sp_mat::const_iterator mu_iter     = mu.begin();
  arma::sp_mat::const_iterator mu_iter_end = mu.end();
  for(; mu_iter != mu_iter_end; ++mu_iter){
    double xval = x(mu_iter.row(),mu_iter.col());
    if(!arma::is_finite(xval)){
      res(mu_iter.row(),mu_iter.col()) = 1 ;
    }
    else {
      res(mu_iter.row(),mu_iter.col()) = R::dnbinom_mu(xval,size_dnb,(*mu_iter),1);
    }
    
  }
  return res;
}
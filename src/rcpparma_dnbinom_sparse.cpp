// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp; 
using namespace arma;

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>

static int NBthreads = -1;

int get_lldist_threads(const int n_profiles) {
  if (NBthreads == -1) {
    // Max allocation of threads equal to 80% of cores
    NBthreads = floor(0.8*omp_get_num_procs());
    
    // Reduce max based on OpenMP settings
    NBthreads = std::min(NBthreads, omp_get_thread_limit());
    NBthreads = std::min(NBthreads, omp_get_max_threads());
  }
  const int ans = n_profiles + 2; // desired number of threads
  return std::min(ans, NBthreads);
}
#endif

//' sum from negative binomial density function
//'
//' Probability density function of the negative binomial distribution (written in C++)
//'
//' @param mat dgCMatrix expression counts
//' @param bgsub vector of background expression per cell
//' @param x numeric expression for reference profiles
//' @param bg numeric background level
//' @param size_dnb int Dispersion parameter
//'
//' @return rowSums for matrix of densities
//' @useDynLib InSituType, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @exportPattern "^[[:alpha:]]+" 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix
  lls_rna(arma::sp_mat& mat, arma::vec& bgsub, arma::mat& x, arma::vec& bg, int& size_dnb) {
    unsigned int K = x.n_cols;
    Rcpp::NumericMatrix res(mat.n_rows, K);
#pragma omp parallel for num_threads(get_lldist_threads(K))
    for (unsigned int k = 0; k < K; k++) {
      const arma::mat::const_col_iterator col_it_begin = x.begin_col(k);
      arma::mat::const_col_iterator col_it = x.begin_col(k);
      const arma::mat::const_col_iterator col_it_end = x.end_col(k);
      const arma::vec s = bgsub / sum(x.col(k));
      for(; col_it != col_it_end; ++col_it) {
        arma::vec::const_iterator s_iter = s.begin();
        arma::vec::const_iterator bg_iter = bg.begin();
        for(; s_iter != s.end(); ++s_iter) {
          double yhat = (*s_iter) * (*col_it) + (*bg_iter);
          int i = s_iter - s.begin();
          int j = col_it - col_it_begin;
          res(i, k) += R::dnbinom_mu(mat(i, j), size_dnb, yhat, 1);
          ++bg_iter;
        }
      }
    }
    return res;
  }

//' sum from Gaussian density function
//'
//' Probability density function of the Gaussian distribution (written in C++)
//'
//' @param mat dgCMatrix expression matrix
//' @param bgsub vector of background expression per cell
//' @param x numeric expression for reference profiles
//' @param xsd numeric expression for reference SD profiles
//' 
//' @return rowSums for matrix of densities
//' @useDynLib InSituType, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @exportPattern "^[[:alpha:]]+" 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix
  lls_protein(arma::mat& mat, arma::vec& bgsub, arma::mat& x, arma::mat& xsd) {
    unsigned int K = x.n_cols;
    Rcpp::NumericMatrix res(mat.n_rows, K);
#pragma omp parallel for num_threads(get_lldist_threads(K))
    for (unsigned int k = 0; k < K; k++) {
      const arma::mat::const_col_iterator col_it_begin = x.begin_col(k);
      arma::mat::const_col_iterator col_it = x.begin_col(k);
      arma::mat::const_col_iterator xsd_iter = xsd.begin_col(k);
      const arma::mat::const_col_iterator col_it_end = x.end_col(k);
      const arma::vec s = bgsub / sum(x.col(k));
      for(; col_it != col_it_end; ++col_it) {
        arma::vec::const_iterator s_iter = s.begin();
        //arma::vec::const_iterator bg_iter = bg.begin();
        for(; s_iter != s.end(); ++s_iter) {
          double yhat = (*s_iter) * (*col_it);
          double sd = (*s_iter) * (*xsd_iter);
          int i = s_iter - s.begin();
          int j = col_it - col_it_begin;
          res(i, k) += R::dnorm(mat(i, j), yhat, sd, 1);
          //++bg_iter;
        }
        ++xsd_iter;
      }
    }
    return res;
  }
 
#include "RcppArmadillo.h"
#include "ct_structure.h"
#include "input_check.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @name ct_structure
//' @title Correlation Thresholded Factor Analysis Structure
//' 
//' @description Generates a factor analysis structure through correlation
//' thresholding.
//'
//' @param R; numeric matrix of correlations.
//' @param tau; numeric scalar of threshold \eqn{\tau} value.
//' @param check; optional logical scalar to indicate if `R` should be
//' checked to be a valid input.
//' 
//' @return `numeric matrix` of the support of the coefficients or support of
//' the factor loadings: \eqn{\mathcal{A}(\Lambda)}.
//' 
//' @references Kim, D. S., & Zhou, Q. (2023). Structure learning of latent
//' factors via clique search on correlation thresholded graphs. *Proceedings
//' of the 40th International Conference on Machine Learning., 202*,
//' 16978–16996. \url{https://proceedings.mlr.press/v202/kim23aa.html}
//' 
//' @export

// [[Rcpp::export]]
arma::umat ct_structure(
  arma::mat const &R, double const &tau, bool const check = true
) {
  
  // I. Preliminaries
  // A. Input Check
  if(check) input_check(R);
  
  // B. Quantities
  const arma::uword p = R.n_rows;
  arma::umat E = arma::abs(R) > tau; // Thresholded Correlation Matrix
  E.diag().ones();
  
  // II. Structure Learning
  std::set<std::vector<int>> lambda_set;
  for(arma::uword i = 0; i < p; i++) {
    
    // A. Independent Maximal Clique Test
    const arma::uvec ne_ind = arma::find(E.col(i) == 1);
    const arma::umat ne_i = E.submat(ne_ind, ne_ind);
    
    // B. Save Structure
    if(arma::accu(ne_i) == ne_i.n_elem && arma::accu(ne_i) > 1) {
      arma::uvec struc_vec(p);
      struc_vec.elem(ne_ind) = arma::uvec(ne_ind.n_elem, arma::fill::ones);
      lambda_set.insert(arma::conv_to<std::vector<int>>::from(struc_vec));
    }
    
  }
  
  // III. Construct Lambda
  arma::umat lambda(p, 0);
  
  // A. Empty Condition
  if(lambda_set.size() == 0) return lambda;
  
  // B. Non-Empty Condition
  for(auto i = lambda_set.begin(); i != lambda_set.end(); i++) {
    lambda.insert_cols(0, arma::conv_to<arma::uvec>::from(*i));
  }
  
  // IV. Return
  return lambda;
  
}

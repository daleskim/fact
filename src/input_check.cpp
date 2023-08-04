#include "RcppArmadillo.h"
#include "ct_structure.h"
// [[Rcpp::depends(RcppArmadillo)]]

void input_check(arma::mat const &R) {
  
  // I. Check Size
  if(R.n_elem == 1) {
    throw std::runtime_error(
      "Correlation matrix size must be greater than 1 x 1."
    );
  }
  
  // II. Check Square
  if(!R.is_square()) {
    throw std::runtime_error("Correlation matrix must be square.");
  }
  
  // III. Check Symmetric
  if(!R.is_symmetric()) {
    throw std::runtime_error("Correlation matrix must be symmetric.");
  }
  
}


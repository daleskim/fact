#include "RcppArmadillo.h"
#include "pare_tau.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @name pare_tau
//' @title Pare Tau
//' 
//' @description Removes redundant \eqn{\tau} values given a correlation matrix
//' \eqn{R}.
//'
//' @param R; numeric matrix of correlations.
//' @param tau; numeric vector of \eqn{\tau} values.
//' 
//' @details Two \eqn{\tau} values are considered redundant if they result in
//' the same thresholded graph. Redundant \eqn{\tau} are removed prior to
//' clique analysis for efficient computation.
//' 
//' @return `numeric vector` of original \eqn{\tau} values with
//' redundant elements removed.
//' 
//' @export

// [[Rcpp::export]]
arma::vec pare_tau(arma::mat const &R, arma::vec const &tau) {
  
  // I. Preliminaries
  const arma::vec tau_arma = tau(arma::find_unique(tau));
  const arma::uword p = R.n_rows;
  const arma::uword pc2 = p * (p - 1) / 2;
  const arma::uword m = tau_arma.n_elem;
  const arma::uvec z = {0};
  
  // II. Create Combined Matrix
  // A. Vector of Values
  const arma::uvec lt_ind = arma::trimatl_ind(arma::size(R), -1);
  const arma::vec cor_vec = arma::abs(R.elem(lt_ind));
  const arma::vec values = arma::join_cols(cor_vec, tau_arma);
  
  // B. Indicator Vectors
  const arma::vec cor_idc(pc2, arma::fill::zeros);
  const arma::vec tau_idc(m, arma::fill::ones);
  const arma::vec idc = arma::join_cols(cor_idc, tau_idc);
  
  // C. Combined Matrix
  arma::mat comb_mat = arma::join_rows(values, idc);
  arma::uvec sort_ind = arma::sort_index(comb_mat.col(0));
  comb_mat = comb_mat.rows(sort_ind);
  
  // III. Pare Tau
  // A. Elements Above Max
  const double R_max = cor_vec.max();
  const arma::uvec above_max_ind = arma::find(
    comb_mat.col(1) == 1 && comb_mat.col(0) >= R_max
  );
  comb_mat.shed_rows(above_max_ind);
  
  // B. Redundant Elements
  arma::uvec tau_ind = arma::find(comb_mat.col(1) == 1);
  arma::uword tau_i = 0;
  
  // 1. Early Return (Single Tau)
  if(tau_ind.n_elem == 1) return comb_mat.submat(tau_ind, z);
  
  // 2. Multiple Tau
  while(true) {
    
    if(tau_ind(tau_i) + 1 == tau_ind(tau_i + 1)) {
      tau_ind.shed_row(tau_i);
    } else {
      tau_i++;
    }
    
    if(tau_i == tau_ind.n_elem - 1) {
      break;
    }
    
  }
  const arma::vec pared_tau = comb_mat.submat(tau_ind, z);
  
  // IV. Return
  return pared_tau;
  
}


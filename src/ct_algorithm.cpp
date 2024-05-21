#include "RcppArmadillo.h"
#include "pare_tau.h"
#include "ct_structure.h"
#include "input_check.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @name ct_algorithm
//' @title Factor Analysis by Correlation Thresholding
//' 
//' @description Generates a set of factor analysis structures by thresholding
//' the correlation graph at a sequence of points.
//'
//' @param R; numeric matrix of correlations.
//' @param tau; optional numeric vector of threshold \eqn{\tau} values.
//' @param k; optional numeric scalar integer indicating minimum clique size
//' to search for. Must be 2 or greater, default is 2.
//' 
//' @details If `tau` is not given, then a set is automatically generated as
//' follows:
//' \itemize{
//'   \item If the number of observed variables is \eqn{\leq 10}, then every
//'   possible threshold is checked. This is at most 45 thresholds.
//'   \item If the number of observed variables is \eqn{> 10}, then 50
//'   evenly spaced thresholds are generated from the interval \eqn{[0, 1]}.
//' }
//' Prior to clique analysis, redundant \eqn{\tau} values are removed for
//' computational efficiency. Therefore, if two \eqn{\tau} values result in
//' identical factor analysis structures, the structure will only be returned
//' once.
//' 
//' @return `list` of `numeric matrix`; each contains the support of the
//' coefficients or support of the factor loadings: \eqn{\mathcal{A}(\Lambda)}.
//' One per (non-redundant) threshold \eqn{\tau}.
//' 
//' @references Kim, D. S., & Zhou, Q. (2023). Structure learning of latent
//' factors via clique search on correlation thresholded graphs. *Proceedings
//' of the 40th International Conference on Machine Learning., 202*,
//' 16978–16996. \url{https://proceedings.mlr.press/v202/kim23aa.html}
//' 
//' @export

// [[Rcpp::export]]
Rcpp::List ct_algorithm(
    const arma::mat &R,
    const Rcpp::Nullable<Rcpp::NumericVector> tau = R_NilValue,
    const int k = 2
) {
  
  // I. Preliminaries
  // A. Check R
  input_check(R);
  
  // B. Check k
  if(k < 2) throw std::runtime_error("k must be 2 or greater.");
  
  // C. Format Tau
  arma::vec tau_arma;
  
  // 1. Tau is Given
  if(tau.isNotNull()) {
    const Rcpp::NumericVector tau_rcpp(tau);
    tau_arma = Rcpp::as<arma::vec>(Rcpp::wrap(tau_rcpp));
    
  // 2. Tau is Automatic (Small p)
  } else if(R.n_rows <= 10 && R.n_rows > 1) {
    
    // a. Preliminaries
    const arma::vec z(1, arma::fill::zeros);
    const arma::uvec lt_ind = arma::trimatl_ind(arma::size(R), -1);
    
    // b. Sort R as Vector
    arma::vec R_vec = R.elem(lt_ind);
    R_vec.insert_rows(0, z);
    R_vec = R_vec(arma::find_unique(R_vec));
    R_vec = arma::sort(R_vec);
    
    // c. Generate Midpoints
    tau_arma.resize(R_vec.n_elem - 1);
    for(arma::uword i = 0; i < tau_arma.n_elem; i++){
      tau_arma(i) = (R_vec(i) + R_vec(i + 1)) / 2;
    }
    
  // 3. Tau is Automatic (Large p)
  } else {
    tau_arma = arma::regspace(0, 0.02, 1);
  }
  
  // B. Pare Tau
  const arma::vec pared_tau = pare_tau(R, tau_arma);
  const arma::uword m = pared_tau.n_elem;
  
  // II. CT Algorithm
  Rcpp::List out_list(0);
  std::set<std::vector<int>> out_set;
  for(arma::uword i = 0; i < m; i++) {
    arma::umat struc = ct_structure(R, pared_tau(i), k, false);
    if(struc.n_cols > 0) {
      const int n_test = out_set.size();
      const std::string tau_name = std::to_string(pared_tau(i));
      arma::umat struc_vec = struc;
      struc_vec.resize(struc.n_rows * struc.n_cols, 1);
      out_set.insert(
        arma::conv_to<std::vector<int>>::from(struc_vec)
      );
      // A. Uniqueness Check
      if(n_test != out_set.size()) out_list[tau_name] = struc;
    }
  }
  
  // III. Return
  return out_list;
  
}

